## Searching the Database

The database is in OPTIMADE format and can be easily searched locally.

**Quick Start:**
1. Install dependencies: "pip install -r dependencies.txt"
2. Install MondoDB Community Server: https://www.mongodb.com/try/download/community
3. Install MongoDB Commandline tools: https://www.mongodb.com/try/download/database-tools
4. Install OPTIMADE into your virtual environment: https://www.optimade.org/optimade-python-tools/latest/INSTALL/
5. Install MongoDB Shell: https://www.mongodb.com/docs/mongodb-shell/install/
6. Activate MongoDB: "sudo systemctl start mongod"
7. Import optimade_result.json into the DB : `mongoimport --db=molecular_dynamics_db --collection=structures --file=optimade_results.json --jsonArray`
8. Start server, referencing your own configuration file: `OPTIMADE_CONFIG_FILE=*absolute path to config mongo_optimade_config.json* python3 -m uvicorn optimade.server.main:app --host 0.0.0.0 --port 5000`
9. Observe the endpoint by searching the following in a web browser: `http://localhost:5000/v1/`

**Bugfix: For if the data displays incorrectly at the endpoint**
10. Run the following CL command: `./flatten.sh`