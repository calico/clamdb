import re
import pandas as pd
import numpy as np
import sys
import os
import feather

def msp_to_df(msp_path, file_encoding = 'UTF-8'):

    '''
    MSP to DataFrame

    Reads a .msp file from msp_path and separates each msp entry into attributes and fragmentation
    data returning two pandas DataFrames - one for attributes, and another for fragmentations.
    '''


    # a global dictionary of entries
    cmpd_attributes = dict()
    cmpd_fragmentation = dict()

    # attributes and fragmentation data for a single entry
    one_cmpd_attributes = dict()

    one_cmpd_frag_ic = dict()
    one_cmpd_frag_label = dict()
    one_cmpd_frag_line = dict()

    # variable to track whether we are in the attributes or fragmentation measurements
    # section of msp file
    tracking_cmpd_attr = True

    cmpd_index = 0
    line_count = 0

    try:
        with open(msp_path, encoding = file_encoding) as msp_file:
            for line in msp_file:
                line_count += 1

                line = line.strip()

                # add attributes & fragments of a compound
                if line:
                    # add attributes
                    if tracking_cmpd_attr is True:
                        split_attr = re.split(':', line, maxsplit = 1)
                        try:
                            one_cmpd_attributes[split_attr[0]] = split_attr[1]
                        except:
                            one_cmpd_attributes[split_attr[0]] = np.nan

                        if split_attr[0] in ["Num Peaks", "Num peaks", "NumPeaks"]:
                            tracking_cmpd_attr = False
                    # add fragments
                    else:
                        split_frag = line.split()

                         # write the frag [mz,ic] pair as a dictionary
                        try:
                            frag_mz = float(split_frag[0])
                        except:
                            frag_mz = np.nan # missing value

                        try:
                            frag_ic = float(split_frag[1])
                        except:
                            frag_ic = np.nan # missing value

                        one_cmpd_frag_ic[frag_mz] = frag_ic

                        # add a label if one is provided
                        if len(split_frag) >= 3:
                            one_cmpd_frag_label[frag_mz] = " ".join(split_frag[2:len(split_frag)])
                        else:
                            one_cmpd_frag_label[frag_mz] = ""

                        # save the frag input line
                        one_cmpd_frag_line[frag_mz] = line

                # compound ends at blank line, save compound attributes and fragments and start new compound
                else:
                    cmpd_attributes[cmpd_index] = pd.DataFrame({'attribute' : list(one_cmpd_attributes.keys()),
                                                                'value' : list(one_cmpd_attributes.values())})
                    cmpd_attributes[cmpd_index]['entry'] = cmpd_index
                    one_cmpd_attributes = dict()

                    cmpd_fragmentation[cmpd_index] = pd.DataFrame({'mz' : list(one_cmpd_frag_ic.keys()),
                                                                   'ic' : list(one_cmpd_frag_ic.values()),
                                                                   'label' : list(one_cmpd_frag_label.values()),
                                                                   'line' : list(one_cmpd_frag_line.values())})
                    cmpd_fragmentation[cmpd_index]['entry'] = cmpd_index
                    one_cmpd_frag_ic = dict()
                    one_cmpd_frag_label = dict()
                    one_cmpd_frag_line = dict()

                    cmpd_index += 1
                    tracking_cmpd_attr = True
                    
                    if cmpd_index % 50000 == 0:
                        print(str(cmpd_index) + ' records parsed')

            cmpd_attributes = pd.concat(cmpd_attributes.values(), ignore_index = True)
            cmpd_fragmentation = pd.concat(cmpd_fragmentation.values(), ignore_index = True)


            return([cmpd_attributes, cmpd_fragmentation])
    except UnicodeDecodeError:
        print("unicode error on line " + str(line_count))

def mgf_to_df(msp_path, file_encoding = 'UTF-8'):

    '''
    MSP to DataFrame

    Reads a .msp file from msp_path and separates each msp entry into attributes and fragmentation
    data returning two pandas DataFrames - one for attributes, and another for fragmentations.
    '''


    # a global dictionary of entries
    cmpd_attributes = dict()
    cmpd_fragmentation = dict()

    # attributes and fragmentation data for a single entry
    one_cmpd_attributes = dict()

    one_cmpd_frag_ic = dict()
    one_cmpd_frag_label = dict()
    one_cmpd_frag_line = dict()

    # variable to track whether we are in the attributes or fragmentation measurements
    # section of msp file
    tracking_cmpd_attr = True

    cmpd_index = 0
    line_count = 0

    try:
        with open(msp_path, encoding = file_encoding) as msp_file:
            for line in msp_file:
                line_count += 1

                line = line.strip()

                # add attributes & fragments of a compound
                if line:
                    if line in ["BEGIN IONS", "END IONS"]:
                        next
                    # add attributes
                    elif tracking_cmpd_attr is True:
                        split_attr = re.split('=', line, maxsplit = 1)
                        try:
                            one_cmpd_attributes[split_attr[0]] = split_attr[1]
                        except:
                            one_cmpd_attributes[split_attr[0]] = np.nan

                        if split_attr[0] in ["SCANS"]:
                            tracking_cmpd_attr = False
                    # add fragments
                    else:
                        split_frag = re.split(r'\t+', line)
                        
                         # write the frag [mz,ic] pair as a dictionary
                        try:
                            frag_mz = float(split_frag[0])
                        except:
                            frag_mz = np.nan # missing value

                        try:
                            frag_ic = float(split_frag[1])
                        except:
                            frag_ic = np.nan # missing value

                        one_cmpd_frag_ic[frag_mz] = frag_ic

                        # add a label if one is provided
                        if len(split_frag) >= 3:
                            one_cmpd_frag_label[frag_mz] = " ".join(split_frag[2:len(split_frag)])
                        else:
                            one_cmpd_frag_label[frag_mz] = ""
                        
                        # save the frag input line
                        one_cmpd_frag_line[frag_mz] = line

                # compound ends at blank line, save compound attributes and fragments and start new compound
                else:
                    cmpd_attributes[cmpd_index] = pd.DataFrame({'attribute' : list(one_cmpd_attributes.keys()),
                                                                'value' : list(one_cmpd_attributes.values())})
                    cmpd_attributes[cmpd_index]['entry'] = cmpd_index
                    one_cmpd_attributes = dict()

                    cmpd_fragmentation[cmpd_index] = pd.DataFrame({'mz' : list(one_cmpd_frag_ic.keys()),
                                                                   'ic' : list(one_cmpd_frag_ic.values()),
                                                                   'label' : list(one_cmpd_frag_label.values()),
                                                                   'line' : list(one_cmpd_frag_line.values())})
                    cmpd_fragmentation[cmpd_index]['entry'] = cmpd_index
                    one_cmpd_frag_ic = dict()
                    one_cmpd_frag_label = dict()
                    one_cmpd_frag_line = dict()

                    cmpd_index += 1
                    tracking_cmpd_attr = True
                    
                    if cmpd_index % 50000 == 0:
                        print(str(cmpd_index) + ' records parsed')

            cmpd_attributes = pd.concat(cmpd_attributes.values(), ignore_index = True)
            cmpd_fragmentation = pd.concat(cmpd_fragmentation.values(), ignore_index = True)


            return([cmpd_attributes, cmpd_fragmentation])
    except UnicodeDecodeError:
        print("unicode error on line " + str(line_count))

def msp_feather_save(msp_path, cmpd_attributes, cmpd_fragmentation):
    save_file_noext = os.path.splitext(msp_path)[0]
    feather.write_dataframe(cmpd_attributes, save_file_noext + '_attr.feather')
    feather.write_dataframe(cmpd_fragmentation, save_file_noext + '_frag.feather')

def format_msp_and_save(msp_path, file_encoding):
    file_extension = os.path.basename(msp_path).split(".")[1]
    if file_extension in ['msp', 'Msp', 'MSP']:
        cmpd_attributes, cmpd_fragmentations = msp_to_df(msp_path, file_encoding)
    elif file_extension in ['mgf', 'Mgf', 'MGF']:
        cmpd_attributes, cmpd_fragmentations = mgf_to_df(msp_path, file_encoding)
    else:
        sys.exit("unknown file extension")
    msp_feather_save(msp_path, cmpd_attributes, cmpd_fragmentations)

if __name__ == "__main__":

    msp_path = sys.argv[1]

    if len(sys.argv) == 3:
      file_encoding = sys.argv[2]
    else:
      file_encoding = "UTF-8"
    print("processing msp or mgf file at " + msp_path + " with " + file_encoding + " encoding\n")

    format_msp_and_save(msp_path, file_encoding)

    print("All Processes Successfully Completed!")
