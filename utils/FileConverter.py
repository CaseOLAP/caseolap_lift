import json

class FileConversion:
    def __init__(self, input_file, output_file):
        self.input_file = input_file
        self.output_file = output_file


    def parse_bool(self, bool_text):
        return bool_text.lower().startswith('y') or bool_text.lower().startswith('t')

    def isfloat(self,num):
        try:
            float(num)
            return True
        except ValueError:
            return False

    def json_to_txt(self):
        # open the json file
        json_file = open(self.input_file, 'r')
        json_data = json.load(json_file)

        # open the txt file
        text_file = open(self.output_file, 'w')

        # write to text file
        for k, v in json_data.items():

            if type(v) == list:
                # for list of lists, handle all lists seperately
                if type(v[0]) == list:
                    for list_idx in range(len(v)):
                        list_text = " ".join(v[list_idx])
                        text = "{key}: {value}\n".format(key = k, value = list_text)
                        text_file.write(text)
                # if it is not a list of lists, just print all contents out normally
                else:
                    # handle lists of type int seperately 
                    if type(v[0]) == int:
                        list_text = " ".join(map(str, v))
                    else:
                        # list_text = " ".join(str(v))
                        list_text = " ".join(v)
                    text = "{key}: {value}\n".format(key = k, value = list_text)
                    text_file.write(text)
            # if the type is not a list, then just print out the values
            else:
                text = "{key}: {value}\n".format(key = k, value = v)
                text_file.write(text)
        
        json_file.close()
        text_file.close()

    def txt_to_json(self):

        # open the txt file
        text = open(self.input_file)
        lines = text.readlines()

        # determine keys
        keys = {}
        for line in lines:
            key = line.split(":")[0]
            if key in keys:
                keys[key] += 1
            else:
                keys[key] = 1

        # save res as dict to be transformed into json
        res = {}

        idx = 0
        for k, v in keys.items():


            # if there is more than one apperance of a key, we need to make a list of lists
            if v > 1:
                items = []
                i = 0
                while i < v:
                    split = lines[idx].split(":")
                    line_text = ":".join(split[1:]).strip()
                    list_text = line_text.split(" ")
                    # list_text = lines[idx].split(":")[1].split(" ")
                    # remove blanks and newlines
                    # list_text = [l.strip("\n") for l in list_text if len(l) > 0]
                    items.append(list_text)
                    i += 1
                    idx += 1
                res[k] = items
            else:
                # if v = 1, then we just need to append whatever is there
                split = lines[idx].split(":")
                line_text = ":".join(split[1:]).strip()
                line_text = line_text.replace("\n", "")

                # determine if k is a bool field
                is_string_list = " " in line_text
                is_bool = 'include' in k
                is_int = line_text.isdigit()
                is_float = self.isfloat(line_text)

                out = line_text

                if is_string_list:
                    out = out.split(" ")
                elif is_bool:
                    out = self.parse_bool(out)
                elif is_int:
                    out = int(out)
                elif is_float:
                    out = float(out)
                elif k == 'cellular_components':
                    if type(out) != list:
                        out = [out]
                elif out == 'None':
                    out = None

                res[k] = out

                idx += 1
        
        # write to files
        # with open(self.output_file, "w") as outfile:
        #     json.dump(res, outfile)

        # close files
        text.close()
        # outfile.close()

        return res

