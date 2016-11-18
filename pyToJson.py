# sauvegarde d'un objet en JSON

with open("obj.json", "w") as outfile:
    json.dump(objpy, outfile)

print("job is saved")

# Lecture du fichier JSON

with open("obj.json", "r") as infile:
    objpy = json.load(infile)

