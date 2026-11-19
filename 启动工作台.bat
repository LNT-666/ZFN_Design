@echo off
chcp 65001 >nul
setlocal
title ZFN 设计工作台 —— 关掉本窗口即停止服务
pushd "%~dp0"

rem 需要一个 Python 3.10 或更新的解释器（代码里用了 3.10 的语法）。
rem 判定方式：import sys 并按 version_info 决定退出码，0 = 够新。

set "PYEXE="
set "PYVER=import sys;raise SystemExit(0 if sys.version_info>=(3,10) else 1)"

rem 1) 环境变量优先（想指定别的 Python 就 set ZFN_PYTHON=...）
if defined ZFN_PYTHON (
  if exist "%ZFN_PYTHON%" (
    "%ZFN_PYTHON%" -c "%PYVER%" >nul 2>nul
    if not errorlevel 1 set "PYEXE=%ZFN_PYTHON%"
  )
)

rem 2) Windows 的 py 启动器
if not defined PYEXE (
  py -3 -c "%PYVER%" >nul 2>nul
  if not errorlevel 1 set "PYEXE=py -3"
)

rem 3) PATH 里的 python
if not defined PYEXE (
  python -c "%PYVER%" >nul 2>nul
  if not errorlevel 1 set "PYEXE=python"
)

rem 4) 常见安装位置（python.org 安装包默认装到这两个地方，但常常没加进 PATH）
if not defined PYEXE (
  for /d %%D in ("%LOCALAPPDATA%\Programs\Python\Python3*") do (
    if not defined PYEXE if exist "%%D\python.exe" (
      "%%D\python.exe" -c "%PYVER%" >nul 2>nul
      if not errorlevel 1 set "PYEXE=%%D\python.exe"
    )
  )
)
if not defined PYEXE (
  for /d %%D in ("C:\Python3*") do (
    if not defined PYEXE if exist "%%D\python.exe" (
      "%%D\python.exe" -c "%PYVER%" >nul 2>nul
      if not errorlevel 1 set "PYEXE=%%D\python.exe"
    )
  )
)

rem 5) 开发机上 DSH 自带的解释器（这台机器打包时用得上）
if not defined PYEXE (
  if exist "%USERPROFILE%\.dsh\dsh-runtimes\dsh-primary-runtime\dependencies\python\python.exe" (
    "%USERPROFILE%\.dsh\dsh-runtimes\dsh-primary-runtime\dependencies\python\python.exe" -c "%PYVER%" >nul 2>nul
    if not errorlevel 1 set "PYEXE=%USERPROFILE%\.dsh\dsh-runtimes\dsh-primary-runtime\dependencies\python\python.exe"
  )
)

if not defined PYEXE (
  echo.
  echo   [x] 没找到 Python 3.10 或更新的版本。
  echo.
  echo       装一个就行：https://www.python.org/downloads/
  echo       安装时勾选 "Add python.exe to PATH"，装完重新双击本文件。
  echo.
  echo       如果已经装了但没被找到，可以设一个环境变量指过去：
  echo         set ZFN_PYTHON=C:\路径\python.exe
  echo.
  pause
  popd
  exit /b 1
)

%PYEXE% -c "import zfn_design" >nul 2>nul
if errorlevel 1 (
  echo.
  echo   [x] 用 %PYEXE% 找不到 zfn_design 模块。
  echo       请确认本文件就在解压出来的 ZFN_Design 目录里（旁边应该能看到
  echo       zfn_design 文件夹），然后重新双击。
  echo.
  pause
  popd
  exit /b 1
)

cls
echo ============================================================
echo    ZFN 设计工作台
echo.
echo    解释器: %PYEXE%
echo    浏览器会自动打开；没弹出来就手动访问下面打印的地址。
echo    用完直接关掉这个黑窗口（或按 Ctrl+C），服务就停了。
echo ============================================================
echo.

%PYEXE% -m zfn_design gui %*

echo.
echo 服务已停止。
pause
popd
