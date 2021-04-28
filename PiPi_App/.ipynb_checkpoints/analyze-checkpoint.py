import sqlite3
import json

conn = sqlite3.connect('results.db')
c = conn.cursor()

results = []

c.execute("select * from results;")

# iterate over database rows
for row in c:
  # extract and save info
  article = json.loads(row[0])
  doc_id = article["abstracts-retrieval-response"]["coredata"]["dc:identifier"]
  subjects = [s["$"] for s in article["abstracts-retreival-resonse"]["coredata"]["subject-areas"]["subject-area"]]
  for s in subjects:
    results.append([doc_id, s])


print(results)
