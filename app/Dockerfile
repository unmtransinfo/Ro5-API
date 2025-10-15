FROM python:3.12-slim

WORKDIR /app

COPY app/requirements.txt ./requirements.txt

RUN pip install --no-cache-dir -r requirements.txt

RUN mkdir -p flasgger_static && \
    cp -r /usr/local/lib/python3.12/site-packages/flasgger/ui3/static/* flasgger_static

COPY . .

ENV APP_PORT=8000

CMD echo "RUNTIME: APP_PORT=${APP_PORT}" && gunicorn --bind "0.0.0.0:${APP_PORT}" --reload app:app