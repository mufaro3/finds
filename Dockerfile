FROM python:3.12-slim

# ---- system dependencies (for docs + builds) ----
RUN apt-get update
RUN apt-get install -y make \
    texlive-latex-base texlive-latex-extra \
    texlive-latex-recommended latexmk
RUN rm -rf /var/lib/apt/lists/*

# ---- working directory = project root ----
WORKDIR /app

# ---- install python dependencies ----
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# ---- copy full project ----
COPY . .

# ---- make src importable (important for Sphinx + runtime) ----
ENV PYTHONPATH="/app/src"

# ---- default runtime command ----
CMD ["python", "src/main.py"]
