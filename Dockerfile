FROM python:3.12-slim

# ---- system dependencies (for docs + builds) ----
RUN apt-get update
RUN apt-get install -y make ffmpeg \
    texlive-latex-base texlive-latex-extra \
    texlive-latex-recommended latexmk
RUN rm -rf /var/lib/apt/lists/*

# ---- working directory = project root ----
WORKDIR /fish

# ---- install python dependencies ----
COPY requirements.txt .
RUN pip install -r requirements.txt

# ---- copy full project ----
COPY . .

# ---- make src importable (important for Sphinx + runtime) ----
ENV PYTHONPATH="/fish"

# ---- default runtime command ----
CMD ["python", "-m", "src.main"]
