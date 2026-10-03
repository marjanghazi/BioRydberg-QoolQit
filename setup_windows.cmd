@echo off
setlocal

echo Creating BioRydberg virtual environment...
python -m venv .venv
call .venv\Scripts\activate

echo Upgrading pip...
python -m pip install --upgrade pip

echo Installing requirements...
pip install -r requirements.txt

echo Registering Jupyter kernel...
python -m ipykernel install --user --name biorrydberg --display-name "BioRydberg"

echo.
echo Setup complete.
echo QoolQit version:
python -c "import qoolqit; print(qoolqit.__version__)"

endlocal
