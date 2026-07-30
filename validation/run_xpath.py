from lxml import etree

tree = etree.parse("xsd/sample.xml")
print(tree.xpath("count(//Plant)"))

# Query 1
total_records = tree.xpath("count(//record)")
print(f"Total Records: {int(total_records)}")

# Query 2
anti_inflammatory = tree.xpath(
    "count(//record[contains(pharmacological_activities_compound,'Anti-inflammatory')])"
)
print(f"Anti-inflammatory Records: {int(anti_inflammatory)}")

# Query 3
imppat = tree.xpath(
    "count(//record[source_database='IMPPAT'])"
)
print(f"IMPPAT Records: {int(imppat)}")

# Query 4
garlic = tree.xpath(
    "//record[plant_species='allium sativum']/phytochemical/text()"
)

print("\nGarlic (Allium sativum) Phytochemicals:")

for compound in garlic:
    print("-", compound)