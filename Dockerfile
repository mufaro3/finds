FROM python:3.12-slim

# ---- system dependencies ----
RUN apt-get update && apt-get install -y \
    make \
    ffmpeg \
    texlive-latex-base \
    texlive-latex-extra \
    texlive-latex-recommended \
    latexmk \
    && rm -rf /var/lib/apt/lists/*

# ---- setup a non-root user ----
RUN useradd -u 1000 -m appuser
USER appuser

# ---- working directory ----
WORKDIR /finds

# ---- install Python deps ----
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# ---- copy project ----
COPY . .

# ---- make package importable ----
ENV PYTHONPATH=/finds

# ---- default command ----
CMD ["python3"]
