
call %~dp0\default-repo-path.bat

if not exist %REPO_DIR%\VERSION (
	echo %REPO_DIR%\VERSION not found, %REPO_DIR% is probably not the root of the EDEN repo !
	exit /b 2
)
set /p VERSION=<"%REPO_DIR%\VERSION"

exit /b 0
