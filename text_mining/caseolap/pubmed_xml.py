import xml.etree.ElementTree as ET

def PubmedXML(file_name):
    with open(file_name,"rb") as f:
        context = ET.iterparse(f, events=("start","end"))
        for event, elem in context:
            yield event,elem
            elem.clear()
            
            
def FullTextFromXML(file_name):
    current_vals = {"pmid":None, "pmcid": None, "text":None}
    xml = PubmedXML(file_name)
    flag = False
    for event, elem in xml:
        #found new article
        if elem.tag == "body":
            print("body tag")
        if elem.tag == "article" and event == "start":
            #yield only if its not the first article tag found
            if flag:
                yield current_vals
            flag = True
            current_vals = {"pmid":None, "pmcid": None, "text":None}
        if elem.tag == "article-id" and event == "start":
            if elem.get('pub-id-type') == 'pmid':
                pmid = elem.text
                current_vals["pmid"] = pmid
            if elem.get('pub-id-type') == 'pmc':
                pmcid = elem.text
                current_vals["pmcid"] = pmcid
        if elem.tag == "body" and event == "start":
            current_vals["text"] = ''.join(elem.itertext())
            
def pmids_from_file(file_name):
    with open(file_name,"r") as f:
        pmids = set(line.strip() for line in f.readlines())
        return pmids