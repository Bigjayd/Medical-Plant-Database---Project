print("Script started")

from lxml import etree

print("Loading XML...")
xml_doc = etree.parse("xsd/sample.xml")

print("Loading XSD...")
xsd_doc = etree.parse("xsd/plant.xsd")

print("Building schema...")
schema = etree.XMLSchema(xsd_doc)

print("Validating...")

if schema.validate(xml_doc):
    print("✅ XML is valid!")
else:
    print("❌ XML is NOT valid!")
    print(schema.error_log)

print("Finished")