FROM python:3.12-slim
ENV PYTHONDONTWRITEBYTECODE=1 PYTHONUNBUFFERED=1
RUN pip install --no-cache-dir numpy==2.5.4
WORKDIR /app/viewer
COPY viewer/app.py viewer/index.html viewer/style.css viewer/plot.js viewer/viewer.js ./
USER 1000:1000
EXPOSE 8767
HEALTHCHECK --interval=30s --timeout=5s --start-period=15s --retries=3 CMD python -c "import urllib.request; urllib.request.urlopen('http://127.0.0.1:8767/api/files', timeout=4).read()"
CMD ["python", "app.py", "--host", "0.0.0.0", "--port", "8767"]
