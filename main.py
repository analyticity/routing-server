from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware

app = FastAPI()
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"]
)

@app.get("/")
async def root():
    return {"message": "Hello World"}

@app.post("/find_route_by_coord")
async def find_route_by_coord():
        
    return {"streets_coord": [],
            "route": list(),
            "src_street": "",
            "dst_street": "",
    }
