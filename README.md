# Ro5 API (Flask + RDKit)

A API backend that computes Lipinski Rule of Five (Ro5) descriptors using RDKit.  
Built with Flask, documented with Swagger (`/apidocs`), and containerized with Docker.

---
## Requirements
- Docker
- Docker Compose

---
## Start the service
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
