
call %~dp0\default-repo-path.bat || goto :error

:: add tooling to PATH
call %~dp0\setpath-i686 || goto :error

:: determine VERSION to tag
call %~dp0\get-version || goto :error
set VERSION_PY=%VERSION%

@echo on
pushd %REPO_DIR%
make clean || goto :error_popd
make -j%NUMBER_OF_PROCESSORS% eden wheel CFLAGS_extra="-static" BUILD=release BUILD_STAMP=%VERSION% WHEEL_PLAT=win32 WHEEL_VERSION=%VERSION_PY% || goto :error_popd

popd

@goto :EOF

:error_popd
@popd
:error
@echo ERROR : Last command failed with error code %errorlevel%. Exiting.
@exit /b %errorlevel%
