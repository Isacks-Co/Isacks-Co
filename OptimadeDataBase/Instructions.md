## Searching the Database

The database is in OPTIMADE format and can be easily searched locally.

**Quick Start:**
1. Install dependencies: "pip install -r dependencies.txt"
2. Install MondoDB Community Server https://www.mongodb.com/try/download/community
3. Install MongoDB Commandline tools https://www.mongodb.com/try/download/database-tools
4. Install OPTIMADE https://www.optimade.org/optimade-python-tools/latest/INSTALL/
5. Activate MongoDB with : "sudo systemctl start mongod"
6. Import optimade_result.json into the DB : `mongoimport --db=molecular_dynamics_db --collection=structures --file=optimade_results.json --jsonArray`
7. Start server: `OPTIMADE_CONFIG_FILE=*absolute path to config mongo_optimade_config.json* python3 -m uvicorn optimade.server.main:app --host 0.0.0.0 --port 5000`
8. Observe the endpoint personally by searching the following in a web browser: `http://localhost:5000/v1/`.