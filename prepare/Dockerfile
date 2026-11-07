FROM python:3.12-slim
ENV PYTHONDONTWRITEBYTECODE=1 PYTHONUNBUFFERED=1
RUN pip install --no-cache-dir numpy==2.5.4 pyarrow==26.0.0
WORKDIR /work
COPY export_parquet.py /tools/export_parquet.py
ENTRYPOINT ["python", "/tools/export_parquet.py"]
