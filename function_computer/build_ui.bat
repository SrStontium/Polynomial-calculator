@echo off
rem ================================================================
rem  Build EasyX GUI version of the polynomial calculator.
rem  Toolchain: Dev-Cpp MinGW64 (TDM-GCC) + EasyX for MinGW
rem ================================================================
set "MINGW=C:\Program Files (x86)\Dev-Cpp\MinGW64"
set "ICONOBJ="

rem --- 1) compile the icon resource -------------------------------
rem     windres fails on paths containing spaces, so use the 8.3
rem     short path for the tool itself.
if not exist app.ico goto skipicon
for %%I in ("%MINGW%\bin\windres.exe") do set "WINDRES=%%~sI"
"%WINDRES%" app.rc -o app_icon.o
if errorlevel 1 goto iconfail
goto skipicon
:iconfail
echo [warn] windres could not rebuild app.rc;
echo        using the existing app_icon.o if present.
:skipicon
if exist app_icon.o set "ICONOBJ=app_icon.o"

rem --- 2) compile and link ----------------------------------------
g++ -std=c++11 -O2 -DUNICODE -D_UNICODE -mwindows poly_ui.cpp %ICONOBJ% -o poly_ui.exe ^
    -I"%MINGW%\include" -L"%MINGW%\lib" ^
    -leasyxw -lgdi32 -luser32 -lole32 -luuid -loleaut32 -limm32 -lwinmm

if errorlevel 1 goto fail
echo.
echo [OK] poly_ui.exe was built with the app icon embedded.
goto end

:fail
echo.
echo [FAILED] Compile error. Check the EasyX path in this file:
echo          %MINGW%

:end
pause
