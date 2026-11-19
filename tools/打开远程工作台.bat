@echo off
chcp 65001 >nul
setlocal
title ZFN 设计工作台（远程 ms01:8797）

set "PY=%USERPROFILE%\.dsh\dsh-runtimes\dsh-primary-runtime\dependencies\python\python.exe"
if not exist "%PY%" set "PY=python"
set "TUNNEL=%USERPROFILE%\.dsh\skills\ms01-server\scripts\ms01_tunnel.py"
set "PORT=8797"

if not exist "%TUNNEL%" (
  echo [错误] 找不到隧道脚本：
  echo        %TUNNEL%
  echo.
  pause
  exit /b 1
)

if "%MS_PASS%"=="" (
  echo 请输入 ms01 的登录密码（只用于本次连接，不会保存）：
  set /p "MS_PASS="
)

echo.
echo 正在建立隧道 127.0.0.1:%PORT%  ^-^->  ms01:%PORT% ...
start "" /b "%PY%" "%TUNNEL%" --remote-port %PORT% --local-port %PORT%
timeout /t 4 /nobreak >nul

start "" "http://127.0.0.1:%PORT%/"
echo.
echo 浏览器已打开 http://127.0.0.1:%PORT%/
echo 本窗口请保持打开：关掉它，网页就断了。
echo.
pause
