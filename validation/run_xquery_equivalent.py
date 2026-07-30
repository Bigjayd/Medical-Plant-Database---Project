from collections import Counter
from lxml import etree

tree = etree.parse("xsd/sample.xml")

plants = tree.xpath("//record/plant_species/text()")

counts = Counter(plants)

print("<results>")

for plant, count in sorted(counts.items()):
    if count > 5:
        print(f"  <summary>")
        print(f"    <Plant>{plant}</Plant>")
        print(f"    <PhytochemicalCount>{count}</PhytochemicalCount>")
        print(f"  </summary>")

print("</results>")