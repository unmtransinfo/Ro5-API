# Ro5 API (Flask + RDKit)

A API backend that computes Lipinski Rule of Five (Ro5) descriptors using RDKit.  
Built with Flask, documented with Swagger (`/apidocs`), and containerized with Docker.

---
## Requirements
- Docker
- Docker Compose
---
## Documentation 
The `/apidocs/` page provides detailed information for all available API endpoints.

- **Production:** https://chiltepin.health.unm.edu/ro5/apidocs/
- **Local:** http://localhost:8000/apidocs/

**UI Repository**: [Ro5-UI](https://github.com/unmtransinfo/Ro5-UI)

---
## Setup (Development)
1. Clone the repository
   ```bash
   git clone https://github.com/unmtransinfo/Ro5-API
   ```
2. cd Ro5-API
3. Start the development environment
    ```bash
    docker compose -f compose-development.yml up --build
    ```

---
## Example request
Ethanol:
```bash
curl -s -X POST http://localhost:8000/ro5 -H "Content-Type: application/json" -d '{"smiles":"CCO","vmax":1}' | jq
```
Aspirin:
```bash
curl -s -X POST http://localhost:8000/ro5 -H "Content-Type: application/json" -d '{"smiles":"CC(=O)OC1=CC=CC=C1C(=O)O","vmax":1}' | jq
```
Caffeine:
```bash
curl -s -X POST http://localhost:8000/ro5 -H "Content-Type: application/json" -d '{"smiles":"CN1C(=O)N(C)C(=O)C(N(C)C=N2)=C12","vmax":1}' | jq
```
