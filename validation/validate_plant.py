from lxml import etree

xml_doc = etree.parse("xsd/plant.xml")
xsd_doc = etree.parse("xsd/plant.xsd")

schema = etree.XMLSchema(xsd_doc)

if schema.validate(xml_doc):
    print("✅ plant.xml is valid!")
else:
    print("❌ Validation failed")
    print(schema.error_log)
    