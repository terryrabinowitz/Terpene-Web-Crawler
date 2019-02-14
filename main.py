#!/usr/bin/python
# -*- coding: iso-8859-15 -*-

import time
import urllib2
from bs4 import BeautifulSoup
import httplib
import matplotlib.pyplot as plt
import numpy as np

np.set_printoptions(precision=2, suppress=True)


def make_relative_percent_table(path):
      time_start = time.time()
      strain_master_load_file = path + 'Analytical360_Absolute.txt'
      strain_master_save_file = path + 'Analytical360_Relative.txt'
      with open(strain_master_save_file, 'wb') as g:
            with open(strain_master_load_file, 'rb') as f:
                  header = f.readline()
                  g.write(header)
                  for line in f.readlines():
                        terpenes_new = []
                        parts = line.split()
                        strain_plus_id = parts[0]
                        g.write(strain_plus_id)
                        g.write("\t")
                        terpenes = parts[1:len(parts)]
                        terpenes = [float(i) for i in terpenes]
                        for i in range(0, len(terpenes) - 1):
                              terpenes_new.append(terpenes[i] / terpenes[len(terpenes) - 1])
                        terpenes_new.append(100.00)
                        terpenes_new = [str(i) for i in terpenes_new]
                        terpenes_new = '\t'.join(terpenes_new) + "\n"
                        g.write(terpenes_new)
      time_end = time.time()
      time_diff = time_end - time_start
      print time_diff


def plot_table(path):
      load_file = path + 'strain_terpene_profiles_ALL_duplicates_removed.txt'
      terpenes = []
      with open(load_file, 'rb') as g:
            header = g.readline()
            header = header.split()
            for terp in header:
                  terpenes.append(terp)
            terpenes.pop(0)
            lines = g.readlines()
      data = []
      for line in lines:
            subdata = []
            line = line.split()
            strain = line[0]
            strain_class = strain.split("___")[0]
            if strain_class == 'skunk':
                  for i in range(1, len(line)):
                        subdata.append(float(line[i]))
                  data.append(subdata)
      data = np.asarray(data)
      print "num samples = ", len(data)

      plt.show(block=True)
      plt.figure()
      x = range(len(terpenes))
      x = [(i + 1) * 2 for i in x]
      x = np.asarray(x)
      medianlineprops = dict(linestyle='-', linewidth=2.5, color='blue')
      meanlineprops = dict(linestyle='-', linewidth=2.5, color='red')
      flierprops = dict(markeredgecolor=None, marker='_', color='black')

      plt.boxplot(data, flierprops=flierprops, medianprops=medianlineprops, meanprops=meanlineprops,
                  positions=x, widths=1, showfliers=True, showmeans=True, meanline=True)
      # plt.boxplot(data, positions=x+.2, widths=.1)
      plt.xticks(x, terpenes, rotation=45)
      plt.show()


def make_table_analytical360(path, counter):
      time_start = time.time()
      amount = 10000
      start = (counter * amount) + 1
      url = 'http://archive.analytical360.com/m/archived/'
      terpene_master_save_file = path + 'Analytical360_Archive_Terpenes_' + str(counter) + ".txt"
      strain_master_save_file = path + 'Analytical360_Archive_Strains_' + str(counter) + ".txt"
      strain_master = {}
      terpene_master = set([])
      for i in range(start, start + amount, 1):
            print i
            url_edit = url + str(i)
            try:
                  page = urllib2.urlopen(url_edit)
                  soup = BeautifulSoup(page, "html.parser")
                  list1 = soup.find_all('h3')
                  strain = str(list1[1])
                  strain = strain.replace('<h3>', '')
                  strain = strain.replace('</h3>', '')
                  strain = strain.replace(' ', '_')
                  strain = strain.lower()
                  strain_flag = '#' + str(i)
                  if strain_flag not in strain and str(i) not in strain:
                        if strain not in strain_master:
                              strain = strain + "___" + str(i)
                              strain_master[strain] = {}
                        else:
                              strain_parts = strain.split("___")
                              if len(strain_parts) > 1:
                                    strain = strain_parts[0] + "___" + str(i)
                              else:
                                    strain = strain + "___" + str(i)
                              strain_master[strain] = {}
                        terpene_flag = False
                        for string in soup.stripped_strings:
                              if 'TERPENE-TOTAL' in string:
                                    break
                              if terpene_flag:
                                    if string.startswith('<'):
                                          pass
                                    else:
                                          parts = string.split()
                                          percent = parts[0].replace('%', '')
                                          if 'mg' in percent:
                                                strain_master.pop(strain, None)
                                                break
                                          terpene = parts[1].lower()
                                          terpene_master.add(terpene)
                                          strain_master[strain][terpene] = percent
                              if 'Terpene Profile' in string:
                                    terpene_flag = True
                        if strain_master[strain] == {}:
                              strain_master.pop(strain, None)
            except:
                  pass
      terpene_master = list(terpene_master)
      terpene_master.sort()
      terpene_string = '\t'.join(terpene_master)
      print
      for key in strain_master:
            print key, strain_master[key]
      print
      print terpene_master
      print
      with open(terpene_master_save_file, 'wb') as f:
            for item in terpene_master:
                  f.write("%s\n" % item)
      with open(strain_master_save_file, 'wb') as f:
            f.write("STRAIN\t%s\n" % terpene_string)
            for strain in sorted(strain_master):
                  f.write("%s\t" % strain)
                  for terpene in terpene_master:
                        if terpene in strain_master[strain]:
                              f.write("%s\t" % strain_master[strain][terpene])
                        else:
                              f.write("0.00\t")
                  f.write("\n")
      time_end = time.time()
      time_diff = time_end - time_start
      print time_diff


def join_tables_analytical360(path, start_counter, max_counter):
      time_start = time.time()
      strain_master_save_file = path + 'strain_terpene_profiles_raw.txt'
      with open(strain_master_save_file, 'wb') as g:
            for i in range(start_counter, (max_counter + 1)):
                  strain_master_load_file = path + 'Analytical360_Archive_Strains_' + str(i) + ".txt"
                  with open(strain_master_load_file, 'rb') as f:
                        line = f.readline()
                        if i == start_counter:
                              print line
                              g.write(line)
                        for line in f.readlines():
                              g.write(line)
      time_end = time.time()
      time_diff = time_end - time_start
      print time_diff


def join_tables_analytical360_100k(path):
      time_start = time.time()
      strain_master_save_file = path + 'strain_terpene_profiles_ALL_raw.txt'
      with open(strain_master_save_file, 'wb') as g:
            index = ['100k', '200k', '300k', '400k', '500k', '600k']
            for i in index:
                  strain_master_load_file = path + 'strain_terpene_profiles_duplicates_removed_' + i + ".txt"
                  with open(strain_master_load_file, 'rb') as f:
                        line = f.readline()
                        if i == '100k':
                              g.write(line)
                        for line in f.readlines():
                              g.write(line)
      time_end = time.time()
      time_diff = time_end - time_start
      print time_diff


def sort_terpene_strains(path):
      time_start = time.time()
      strain_master_load_file = path + 'strain_terpene_profiles_raw.txt'
      strain_master_save_file = path + 'strain_terpene_profiles_sorted.txt'
      with open(strain_master_load_file, 'rb') as g:
            line = g.readline()
            lines = g.readlines()
            lines.sort()
      with open(strain_master_save_file, 'wb') as g:
            g.write(line)
            g.writelines(lines)
      time_end = time.time()
      time_diff = time_end - time_start
      print time_diff


def sort_terpene_strains_100k(path):
      time_start = time.time()
      strain_master_load_file = path + 'strain_terpene_profiles_ALL_raw.txt'
      strain_master_save_file = path + 'strain_terpene_profiles_ALL_sorted.txt'
      with open(strain_master_load_file, 'rb') as g:
            line = g.readline()
            lines = g.readlines()
            lines.sort()
      with open(strain_master_save_file, 'wb') as g:
            g.write(line)
            g.writelines(lines)
      time_end = time.time()
      time_diff = time_end - time_start
      print time_diff


def edit_terpene_strains(path):
      time_start = time.time()
      strain_master_load_file = path + 'strain_terpene_profiles_sorted.txt'
      strain_master_save_file = path + 'strain_terpene_profiles_edited.txt'
      bad_list = ['wax', 'shatter', 'hash', 'kief', 'oil', 'bho', 'co2', 'honeycomb', 'solvent', 'sugar', 'blend',
                  'concentrate', 'extract', 'edible', 'butane', 'ethanol', 'decarb,' 'blend', 'trim', 'resin', 'mix',
                  'rso', 'bubble', 'melt', 'alcohol', 'butter', 'budder', 'rosin', 'precursor', 'glass', 'vape', 'rso',
                  'decarb']
      with open(strain_master_save_file, 'wb') as g:
            with open(strain_master_load_file, 'rb') as f:
                  line = f.readline()
                  parts_header = line.split()
                  for item in parts_header:
                        item = item.strip()
                        if item == 'alpha' or item == 'beta':
                              pass
                        else:
                              g.write("%s\t" % item)
                  g.write("TOTAL\n")
                  for line in f.readlines():
                        parts = line.split()
                        strain = parts[0]
                        id = (strain.split("___"))[1]
                        strain = (strain.split("___"))[0]
                        strain = (strain.split("_("))[0]
                        strain = strain.replace("#", "")
                        flag = True
                        for i in bad_list:
                              if i in strain:
                                    flag = False
                        if flag:
                              strain = strain + "___" + id
                              alpha = parts[1]
                              alpha_pinene = parts[2]
                              if alpha != '0':
                                    alpha_pinene = alpha
                              beta = parts[3]
                              beta_pinene = parts[4]
                              if beta != '0':
                                    beta_pinene = beta
                              caryophyllene = parts[5]
                              humulene = parts[6]
                              limonene = parts[7]
                              linalool = parts[8]
                              myrcene = parts[9]
                              ocimene = parts[10]
                              terpinolene = parts[11]
                              total = float(alpha_pinene) + float(beta_pinene) + float(caryophyllene) + float(
                                    humulene) + float(limonene) + float(linalool) + float(myrcene) + float(
                                    ocimene) + float(terpinolene)
                              g.write("%s\t" % strain)
                              g.write("%s\t" % alpha_pinene)
                              g.write("%s\t" % beta_pinene)
                              g.write("%s\t" % caryophyllene)
                              g.write("%s\t" % humulene)
                              g.write("%s\t" % limonene)
                              g.write("%s\t" % linalool)
                              g.write("%s\t" % myrcene)
                              g.write("%s\t" % ocimene)
                              g.write("%s\t" % terpinolene)
                              g.write("%s\n" % str(total))
      time_end = time.time()
      time_diff = time_end - time_start
      print time_diff


def edit_terpene_strains_remove_duplicates(path):
      time_start = time.time()
      strain_master_load_file = path + 'strain_terpene_profiles_edited.txt'
      strain_master_save_file = path + 'strain_terpene_profiles_duplicates_removed.txt'
      strain_master = []
      with open(strain_master_load_file, 'rb') as f:
            header = f.readline()
            for line in f.readlines():
                  parts = line.split()
                  strain_plus_id = parts[0]
                  strain_plus_id_parts = strain_plus_id.split("___")
                  strain = strain_plus_id_parts[0]
                  id = strain_plus_id_parts[1]
                  terpenes = parts[1:len(parts)]
                  strain_master.append([strain, id, terpenes])
      with open(strain_master_save_file, 'wb') as f:
            f.write(header)
            strain_0 = strain_master[0][0] + strain_master[0][1]
            f.write(strain_0)
            for j in strain_master[0][2]:
                  f.write("\t%s" % j)
            f.write("\n")
            for i in range(1, len(strain_master)):
                  strain_current = strain_master[i][0] + "___" + strain_master[i][1]
                  terpenes_current = strain_master[i][2]
                  print_flag = True
                  for j in range(i - 1, -1, -1):
                        terpenes_compare = strain_master[j][2]
                        if terpenes_current == terpenes_compare:
                              print_flag = False
                  if print_flag:
                        f.write(strain_current)
                        for j in strain_master[i][2]:
                              f.write("\t%s" % j)
                        f.write("\n")
      time_end = time.time()
      time_diff = time_end - time_start
      print time_diff


def edit_terpene_strains_remove_duplicates_100k(path):
      time_start = time.time()
      strain_master_load_file = path + 'strain_terpene_profiles_ALL_sorted.txt'
      strain_master_save_file = path + 'strain_terpene_profiles_ALL_duplicates_removed.txt'
      strain_master = []
      with open(strain_master_load_file, 'rb') as f:
            header = f.readline()
            for line in f.readlines():
                  parts = line.split()
                  print parts
                  strain_plus_id = parts[0]
                  strain_plus_id_parts = strain_plus_id.split("___")
                  strain = strain_plus_id_parts[0]
                  id = strain_plus_id_parts[1]
                  terpenes = parts[1:len(parts)]
                  strain_master.append([strain, id, terpenes])
      with open(strain_master_save_file, 'wb') as f:
            f.write(header)
            strain_0 = strain_master[0][0] + strain_master[0][1]
            f.write(strain_0)
            for j in strain_master[0][2]:
                  f.write("\t%s" % j)
            f.write("\n")
            for i in range(1, len(strain_master)):
                  strain_current = strain_master[i][0] + "___" + strain_master[i][1]
                  terpenes_current = strain_master[i][2]
                  print_flag = True
                  for j in range(i - 1, -1, -1):
                        terpenes_compare = strain_master[j][2]
                        if terpenes_current == terpenes_compare:
                              print_flag = False
                  if print_flag:
                        f.write(strain_current)
                        for j in strain_master[i][2]:
                              f.write("\t%s" % j)
                        f.write("\n")
      time_end = time.time()
      time_diff = time_end - time_start
      print time_diff


def make_statistics(path, abs):
      time_start = time.time()
      if abs:
            strain_master_load_file = path + 'Analytical360_Absolute.txt'
            strain_master_save_file_mean = path + 'Analytical360_Absolute_Mean.txt'
            strain_master_save_file_median = path + 'Analytical360_Absolute_Median.txt'
      else:
            strain_master_load_file = path + 'Analytical360_Relative.txt'
            strain_master_save_file_mean = path + 'Analytical360_Relative_Mean.txt'
            strain_master_save_file_median = path + 'Analytical360_Relative_Median.txt'
      strain_master = {}
      with open(strain_master_load_file, 'rb') as f:
            header = f.readline()
            for line in f.readlines():
                  parts = line.split()
                  strain_plus_id = parts[0]
                  strain_plus_id_parts = strain_plus_id.split("___")
                  strain = strain_plus_id_parts[0]
                  terpenes = parts[1:len(parts)]
                  terpenes = np.asarray([float(i) for i in terpenes])
                  if strain in strain_master:
                        strain_master[strain].append(terpenes)
                  else:
                        strain_master[strain] = [terpenes]
            with open(strain_master_save_file_mean, 'wb') as f:
                  headerx = header.split("STRAIN")
                  headerx = "Strain\tNum-Samples" + headerx[1]
                  f.write(headerx)
                  for strain in sorted(strain_master):
                        num_samples = len(strain_master[strain])
                        mean_samples = np.mean(strain_master[strain], axis=0)
                        f.write(strain)
                        f.write("\t")
                        f.write(str(num_samples))
                        for j in mean_samples:
                              f.write("\t%s" % j)
                        f.write("\n")
            with open(strain_master_save_file_median, 'wb') as f:
                  headerx = header.split("STRAIN")
                  headerx = "Strain\tNum-Samples" + headerx[1]
                  f.write(headerx)
                  for strain in sorted(strain_master):
                        num_samples = len(strain_master[strain])
                        median_samples = np.median(strain_master[strain], axis=0)
                        f.write(strain)
                        f.write("\t")
                        f.write(str(num_samples))
                        for j in median_samples:
                              f.write("\t%s" % j)
                        f.write("\n")
      time_end = time.time()
      time_diff = time_end - time_start
      print time_diff


def make_reference_strain_list(path):
      time_start = time.time()
      strain_master_load_file = path + "strain_terpene_profiles_edited.txt"
      strain_master_save_file = path + 'Analytical360_reference_strain_list.txt'
      master_strain_list = []
      with open(strain_master_load_file, 'rb') as f:
            f.readline()
            for line in f.readlines():
                  parts = line.split()
                  master_strain_list.append(parts[0])
      with open(strain_master_save_file, 'wb') as g:
            for item in master_strain_list:
                  g.write("%s\n" % item)
      time_end = time.time()
      time_diff = time_end - time_start
      print time_diff


def make_table_mcr_labs(path):
      time_start = time.time()
      strain_master_file = path + 'MCR_Labs_Strains.txt'
      strain_master_load_file = path + 'MCR_Labs_Manual_Strain_List.txt'
      strain_list_reference_list = set([])
      strain_master = {}
      with open(strain_master_load_file, 'rb') as f:
            for line in f.readlines():
                  strain = line.split("___")[0]
                  strain = strain.replace("_", "-")
                  strain_list_reference_list.add(strain)
      strain_list_reference_list = list(strain_list_reference_list)
      strain_list_reference_list.sort()
      for strain_reference in strain_list_reference_list:
            for i in range(1, 31):
                  strain_reference = strain_reference.strip()
                  if i == 1:
                        url_edit = "/" + strain_reference + "/"
                  else:
                        url_edit = "/" + strain_reference + "-" + str(i) + "/"
                  try:
                        c = httplib.HTTPSConnection('mcrlabs.com')
                        c.request('GET', url_edit)
                        response = c.getresponse()
                        check = response.status
                        if check == 200:
                              url_edit = url_edit.strip('/')
                              strain_master[url_edit] = []
                              data = response.read()
                              flag=False
                              print_flag = False
                              data = data.split("\n")
                              for line in data:
                                    if ('Concentrate' or 'Extract' or 'concentrate' or 'extract' or 'MIP' or 'mip') in line:
                                          strain_master.pop(url_edit)
                                          break
                                    if "<div class=" in line:
                                          flag=False
                                    if flag == True and "<td>" in line:
                                          terpene = line.strip("<td>").strip("</td>").strip("</td> ").lower()
                                          terpene = terpene.replace("\xce\xb1","alpha")
                                          terpene = terpene.replace("\xce\xb2", "beta")
                                          terpene = terpene.replace("\xce\xb3", "gamma")
                                          terpene = terpene.replace("\xce\xb4", "delta")
                                          terpene = terpene.replace("\xcf\x81", "para")
                                    if flag == True and "<span class=" in line:
                                          if 'Not tested' in line:
                                                strain_master.pop(url_edit)
                                                break
                                          else:
                                                if 'Not detected' in line:
                                                      percent = 0.0
                                                else:
                                                      line = line.split("%</span>")[0]
                                                      line = line.split(">")[1]
                                                      percent = float(line)
                                                strain_master[url_edit].append([terpene, percent])
                                                print_flag=True
                                    if "<th>Terpene</th>" in line:
                                          flag=True
                              if print_flag:
                                    print  url_edit, "\t", strain_master[url_edit]
                  except:
                        pass
      with open(strain_master_file, 'wb') as f:
            flag = True
            for strain in sorted(strain_master):
                  if flag:
                        f.write("STRAIN\t")
                        for terp in strain_master[strain]:
                              f.write("%s\t" % terp[0])
                        f.write("\n")
                        flag = False
                  f.write("%s\t" % strain)
                  for terp in strain_master[strain]:
                        f.write("%s\t" % terp[1])
                  f.write("\n")
      time_end = time.time()
      time_diff = time_end - time_start
      print time_diff


def main():
      path = '/Users/terryrabinowitz/Desktop/Cannabis/terpene_web_data/'
      counter = 61
      # make_table_analytical360(path, counter)
      start_counter = 50
      max_counter = 59
      # join_tables_analytical360(path, start_counter, max_counter)
      # sort_terpene_strains(path)
      # edit_terpene_strains(path)
      # edit_terpene_strains_remove_duplicates(path)
      # join_tables_analytical360_100k(path)
      # sort_terpene_strains_100k(path)
      # edit_terpene_strains_remove_duplicates_100k(path)
      #make_relative_percent_table(path)
      #make_statistics(path, abs=False)
      #make_statistics(path, abs=True)

      # make_reference_strain_list(path)
      #make_table_mcr_labs(path)
      plot_table(path)


main()
