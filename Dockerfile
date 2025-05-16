FROM python:3.12

WORKDIR /app

COPY . /app

RUN pip install -r requirements.txt

CMD [ "uvicorn", "main:app", "--port=8001", "--host=0.0.0.0", "--workers=10"]