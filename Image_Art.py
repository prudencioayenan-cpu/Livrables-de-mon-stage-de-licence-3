from flask import Flask, render_template, request, jsonify, send_from_directory, send_file
import os
import io
import base64
import uuid
import numpy as np
from PIL import Image, ImageFilter
import cv2

app = Flask(__name__)

# ── Sécurité de base ──────────────────────────────────────────────────────────
app.config["MAX_CONTENT_LENGTH"] = 25 * 1024 * 1024  # 25 Mo max par fichier envoyé
# ─────────────────────────────────────────────────────────────────────────────

BASE_DIR         = os.path.dirname(os.path.abspath(__file__))
IMAGES_DIR       = os.path.join(BASE_DIR, "images")
UPLOADS_TEMP_DIR = os.path.join(BASE_DIR, "uploads_temp")

ALLOWED_EXT = {".jpg", ".jpeg", ".png", ".tif", ".tiff", ".bmp", ".gif", ".webp"}


def get_images():
    files = []
    for f in sorted(os.listdir(IMAGES_DIR)):
        if os.path.splitext(f)[1].lower() in ALLOWED_EXT:
            files.append(f)
    return files


def est_une_vraie_image(filepath):
    """Vérifie que le contenu du fichier est bien une image valide (pas juste l'extension)."""
    try:
        with Image.open(filepath) as img:
            img.verify()
        return True
    except Exception:
        return False


def open_image_rgb(path):
    img = Image.open(path)
    if hasattr(img, 'n_frames') and img.n_frames > 1:
        img.seek(0)
    if img.mode in ('I;16', 'I'):
        arr = np.array(img, dtype=np.float32)
        arr = ((arr - arr.min()) / (arr.max() - arr.min() + 1e-8) * 255).astype(np.uint8)
        img = Image.fromarray(arr)
    return img.convert("RGB")


def image_vers_jpeg_bytes(img_pil):
    buf = io.BytesIO()
    img_pil.save(buf, format="JPEG", quality=90)
    buf.seek(0)
    return buf


# ── VERSION 1 : Maël Zami ────────────────────────────────────────────────────
# Flou gaussien + quantification couleurs + Find Edges
def vitrail_v1_mael(img_pil, radius, num_colors, lead_thick):
    arr = np.array(img_pil, dtype=np.uint8)
    # Enhance contrast
    for c in range(3):
        ch = arr[:, :, c]
        p_low, p_high = np.percentile(ch, 0.5), np.percentile(ch, 99.5)
        arr[:, :, c] = np.clip((ch - p_low) / (p_high - p_low + 1e-8) * 255, 0, 255).astype(np.uint8)
    # Flou gaussien
    blurred = cv2.GaussianBlur(arr, (0, 0), sigmaX=max(1, radius))
    # Quantification couleurs
    img_q = Image.fromarray(blurred).quantize(colors=num_colors, method=Image.Quantize.MEDIANCUT).convert("RGB")
    arr_q = np.array(img_q, dtype=np.uint8)
    # Contours
    gray  = cv2.cvtColor(arr_q, cv2.COLOR_RGB2GRAY)
    edges = cv2.Canny(gray, 20, 80)
    if lead_thick > 1:
        k = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (lead_thick*2-1, lead_thick*2-1))
        edges = cv2.dilate(edges, k)
    mask = cv2.bitwise_not(edges)
    mask_rgb = cv2.cvtColor(mask, cv2.COLOR_GRAY2RGB)
    return Image.fromarray(np.minimum(arr_q, mask_rgb).astype(np.uint8))


# ── VERSION 3 : Cléo Thury ───────────────────────────────────────────────────
# Find Maxima + couleur moyenne + boost saturation HSV + choix teinte
def vitrail_v3_cleo(img_pil, sensibilite, precision, intensite, teinte):
    arr = np.array(img_pil, dtype=np.uint8)
    # Style teinte
    if teinte == "Niveaux de gris":
        gray3 = cv2.cvtColor(arr, cv2.COLOR_RGB2GRAY)
        arr = cv2.cvtColor(gray3, cv2.COLOR_GRAY2RGB)
    # Boost saturation HSV (équivalent HSB Stack * intensite)
    if intensite != 1.0:
        hsv = cv2.cvtColor(arr, cv2.COLOR_RGB2HSV).astype(np.float32)
        hsv[:, :, 1] = np.clip(hsv[:, :, 1] * intensite, 0, 255)
        arr = cv2.cvtColor(hsv.astype(np.uint8), cv2.COLOR_HSV2RGB)
    # Segmentation superpixels
    region_size = max(5, int((101 - sensibilite) / 2))
    slic = cv2.ximgproc.createSuperpixelSLIC(arr, region_size=region_size, ruler=float(precision))
    slic.iterate(10)
    labels = slic.getLabels()
    result = np.zeros_like(arr)
    for lbl in range(labels.max() + 1):
        mask = (labels == lbl)
        for c in range(3):
            result[:, :, c][mask] = int(arr[:, :, c][mask].mean())
    # Plomb
    contour_mask = slic.getLabelContourMask()
    contour_inv  = cv2.bitwise_not(contour_mask)
    contour_rgb  = cv2.cvtColor(contour_inv, cv2.COLOR_GRAY2RGB)
    return Image.fromarray(np.minimum(result, contour_rgb).astype(np.uint8))


# ── Mosaïque ─────────────────────────────────────────────────────────────────
def traiter_mosaique(img_pil, dim, bord, brill):
    arr = np.array(img_pil, dtype=np.uint8)
    h, w = arr.shape[:2]
    result = np.zeros((h, w, 3), dtype=np.uint8)
    cols, rows = w // dim, h // dim
    offset, size = bord // 2, dim - bord
    if size <= 0:
        return img_pil
    for y in range(rows):
        for x in range(cols):
            px, py = x * dim, y * dim
            cx, cy = min(px + dim // 2, w - 1), min(py + dim // 2, h - 1)
            r, g, b = arr[cy, cx]
            r = int(min(255, r * brill))
            g = int(min(255, g * brill))
            b = int(min(255, b * brill))
            x1, y1 = px + offset, py + offset
            result[y1:y1+size, x1:x1+size] = [r, g, b]
    return Image.fromarray(result)


# ── Routes Flask ──────────────────────────────────────────────────────────────
@app.route("/")
def index():
    return render_template("index.html", images=get_images())


@app.route("/images/<path:filename>")
def serve_image(filename):
    ext = os.path.splitext(filename)[1].lower()
    if ext in {".tif", ".tiff"}:
        try:
            img = open_image_rgb(os.path.join(IMAGES_DIR, filename))
            return send_file(image_vers_jpeg_bytes(img), mimetype="image/jpeg")
        except Exception as e:
            return f"Erreur : {e}", 500
    return send_from_directory(IMAGES_DIR, filename)


@app.route("/uploads_temp/<path:filename>")
def serve_uploaded_image(filename):
    """Sert la vignette d'une image temporairement téléversée par l'utilisateur."""
    ext = os.path.splitext(filename)[1].lower()
    filepath = os.path.join(UPLOADS_TEMP_DIR, filename)
    if not os.path.isfile(filepath):
        return "Fichier introuvable.", 404
    if ext in {".tif", ".tiff"}:
        try:
            img = open_image_rgb(filepath)
            return send_file(image_vers_jpeg_bytes(img), mimetype="image/jpeg")
        except Exception as e:
            return f"Erreur : {e}", 500
    return send_from_directory(UPLOADS_TEMP_DIR, filename)


@app.route("/upload", methods=["POST"])
def upload():
    """Reçoit une image envoyée par l'utilisateur et la stocke temporairement."""
    if "fichier" not in request.files:
        return jsonify({"ok": False, "erreur": "Aucun fichier reçu."})

    fichier = request.files["fichier"]
    if fichier.filename == "":
        return jsonify({"ok": False, "erreur": "Nom de fichier vide."})

    ext = os.path.splitext(fichier.filename)[1].lower()
    if ext not in ALLOWED_EXT:
        return jsonify({"ok": False, "erreur": f"Format non supporté ({ext}). Formats acceptés : {', '.join(sorted(ALLOWED_EXT))}."})

    os.makedirs(UPLOADS_TEMP_DIR, exist_ok=True)

    # Nom unique pour éviter les collisions entre utilisateurs
    nom_unique = f"{uuid.uuid4().hex}{ext}"
    chemin_temp = os.path.join(UPLOADS_TEMP_DIR, nom_unique)
    fichier.save(chemin_temp)

    # Vérifie que c'est une vraie image et pas un fichier malveillant renommé
    if not est_une_vraie_image(chemin_temp):
        os.remove(chemin_temp)
        return jsonify({"ok": False, "erreur": "Le fichier envoyé n'est pas une image valide."})

    return jsonify({"ok": True, "nom": nom_unique, "nom_original": fichier.filename})


@app.route("/lancer", methods=["POST"])
def lancer():
    data       = request.json
    image_name = data.get("image")
    programme  = data.get("programme")
    params     = data.get("params", {})
    est_upload = data.get("estUpload", False)

    # Cherche l'image dans le bon dossier selon son origine
    dossier_source = UPLOADS_TEMP_DIR if est_upload else IMAGES_DIR
    image_path = os.path.join(dossier_source, image_name)
    if not os.path.isfile(image_path):
        return jsonify({"ok": False, "erreur": "Image introuvable."})

    try:
        img = open_image_rgb(image_path)

        if programme == "mosaique":
            resultat = traiter_mosaique(img,
                int(params.get("dim", 20)),
                int(params.get("bord", 2)),
                float(params.get("brill", 1.5)))

        elif programme == "vitrail_v1":
            resultat = vitrail_v1_mael(img,
                int(params.get("radius", 4)),
                int(params.get("numColors", 12)),
                int(params.get("leadThick", 2)))

        elif programme == "vitrail_v3":
            resultat = vitrail_v3_cleo(img,
                int(params.get("sensibilite", 50)),
                int(params.get("precision", 3)),
                float(params.get("intensite", 1.0)),
                params.get("teinte", "Image originale"))

        else:
            return jsonify({"ok": False, "erreur": "Programme inconnu."})

        # Génération du TIF entièrement en mémoire (rien n'est écrit sur le disque)
        buf_tif = io.BytesIO()
        resultat.save(buf_tif, format="TIFF")
        tif_b64 = base64.b64encode(buf_tif.getvalue()).decode("utf-8")

        # Aperçu JPEG pour affichage immédiat dans l'interface
        buf_jpeg = image_vers_jpeg_bytes(resultat)
        apercu   = base64.b64encode(buf_jpeg.getvalue()).decode("utf-8")

        base_name   = os.path.splitext(image_name)[0]
        output_name = f"{programme}_{base_name}.tif"

        return jsonify({
            "ok": True,
            "message": "Traitement terminé !",
            "fichier": output_name,
            "apercu": apercu,
            "tif_b64": tif_b64
        })

    except Exception as e:
        return jsonify({"ok": False, "erreur": str(e)})


@app.route("/telecharger", methods=["POST"])
def telecharger():
    """Reçoit le TIF en base64 (déjà généré par /lancer) et le propose en téléchargement, sans rien stocker sur le serveur."""
    data        = request.json
    tif_b64     = data.get("tif_b64")
    nom_fichier = data.get("fichier", "resultat.tif")

    if not tif_b64:
        return "Aucune donnée à télécharger.", 400

    buf = io.BytesIO(base64.b64decode(tif_b64))
    buf.seek(0)
    return send_file(buf, mimetype="image/tiff", as_attachment=True, download_name=nom_fichier)


if __name__ == "__main__":
    os.makedirs(IMAGES_DIR, exist_ok=True)
    os.makedirs(UPLOADS_TEMP_DIR, exist_ok=True)
    #app.run(debug=False, port=5000)
    app.run(host="0.0.0.0", port=5000, ssl_context="adhoc")