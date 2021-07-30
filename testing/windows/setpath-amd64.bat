@echo off
call %~dp0\default-tooling-path.bat

call :add_to_path %TOOLING_DIR%\coreutils\bin
call :add_to_path %TOOLING_DIR%\flexbison
call :add_to_path %TOOLING_DIR%\make
call :add_to_path %TOOLING_DIR%\xxd

call :add_to_path %TOOLING_DIR%\mingw64\bin
call :add_to_path %TOOLING_DIR%\python64

@exit /b 0

:add_to_path
:: %~f1 expands %1 to a fully qualified path name
:: echo %~f1
:: Check if this is already added to PATH,, since variables persist in this shell
@echo ";%path%;" | %SYSTEMROOT%\System32\find /C /I ";%~f1;" > nul && goto :eof
REM path %~f1;%path%
if %errorlevel%==1 path %~f1;%path%
goto :eof
