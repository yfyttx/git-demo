@echo off
setlocal

powershell -NoProfile -ExecutionPolicy Bypass -File "%~dp0generate_matrices.ps1" %*

echo.
pause
