FROM python:3.12

WORKDIR /app

COPY . /app

RUN pip install --no-cache-dir -r requirements.txt

EXPOSE 8001

CMD [ "uvicorn", "main:app", "--port=8001", "--host=0.0.0.0", "--workers=1"]