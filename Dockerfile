# --- Base με micromamba (conda-forge) ---
FROM mambaorg/micromamba:1.5.8-focal

# Ενεργοποίηση conda env σε ΚΑΘΕ RUN
ARG MAMBA_DOCKERFILE_ACTIVATE=1
SHELL ["/usr/local/bin/_dockerfile_shell.sh"]

ENV DEBIAN_FRONTEND=noninteractive \
    PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    MAMBA_NO_SHORT_ERRORS=1

# 1) Conda deps (ΧΩΡΙΣ pip μέσα στο environment.yml)
COPY environment.yml /tmp/environment.yml
RUN micromamba install -y -n base -f /tmp/environment.yml && \
    micromamba clean --all --yes

# 2) Pip deps σε ξεχωριστό βήμα (πιο σταθερό / καθαρά errors)
COPY requirements-pip.txt /tmp/requirements-pip.txt
RUN pip install --no-cache-dir -r /tmp/requirements-pip.txt

# 3) Κώδικας εφαρμογής
WORKDIR /app
COPY . /app

# Streamlit
ENV STREAMLIT_SERVER_HEADLESS=true \
    STREAMLIT_SERVER_PORT=8501 \
    STREAMLIT_BROWSER_GATHER_USAGE_STATS=false

EXPOSE 8501
# ΠΡΙΝ (μην το κρατήσεις)
# CMD ["streamlit","run","app.py","--server.address=0.0.0.0","--server.port=${STREAMLIT_SERVER_PORT}"]

# ΜΕΤΑ (σταθερά 8501 μέσα στο container)
CMD ["streamlit","run","app.py","--server.address=0.0.0.0","--server.port=8501"]


