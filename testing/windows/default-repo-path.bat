
@if "%REPO_DIR%"=="" (
	set REPO_DIR=%~dp0\..\..\
	:: make content or printing of REPO_DIR= more elegant if you can i dare you
	@echo:REPO_DIR was not defined, set to %~dp0\..\..\ )

exit /b 0
