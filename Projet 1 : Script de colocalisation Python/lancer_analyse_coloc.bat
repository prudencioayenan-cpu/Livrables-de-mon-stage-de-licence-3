@echo off
title Analyse STORM - ImageUP
echo Lancement du logiciel d'analyse STORM...

:: Force Windows a se placer dans le dossier actuel ou se trouve ce fichier .bat
cd /d "%~dp0"

:: Lance le script d'analyse Python local
start python analyse_colocalisation_STORM.py

echo.
echo Le script a ete lance avec succes.
pause