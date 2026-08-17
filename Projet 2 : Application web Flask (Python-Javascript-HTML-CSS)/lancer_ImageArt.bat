
@echo off
:: Se déplace dans le dossier où est stocké ce fichier .bat
cd /d "%~dp0"
:: Lance le script Python local
start python Image_Art.py
timeout /t 3 /nobreak
:: Ouvre le navigateur local
start https://127.0.0.1:5000