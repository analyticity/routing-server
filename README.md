# Routing Server

Routing server integrating OpenStreetMap map data with Waze historical traffic information.

## Running the server

### Prerequisites

- Database credentials in `.env` file
  - `DB_HOST`
  - `DB_PORT`
  - `DB_USER`
  - `DB_PASSWORD`
  - `DB_NAME`
- Python 3.12+ for local setup
- Docker and Docker Compose for Docker setup

### Local setup
```
pip install -r requirements.txt
chmod +x run.sh
./run.sh
```
   
### Docker setup
```
docker compose up -d
```
