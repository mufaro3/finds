# Use Python 3 base image
FROM python:3.12-slim

# Set working directory inside container
WORKDIR /src

# Copy dependency list first (better caching)
COPY requirements.txt .

# Install dependencies
RUN pip3 install --no-cache-dir -r requirements.txt

# Copy project files
COPY . .

# Run application
CMD ["python", "src/main.py"]