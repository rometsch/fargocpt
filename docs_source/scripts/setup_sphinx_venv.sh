if [[ ! -d .venv ]]; then
    echo "Setting up new python virtual env"
    python3 -m venv .venv
    source .venv/bin/activate
    pip install --upgrade pip
    pip install --upgrade -r requirements.txt
fi

source .venv/bin/activate