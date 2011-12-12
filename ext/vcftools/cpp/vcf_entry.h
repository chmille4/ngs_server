/*
 * vcf_entry.h
 *
 *  Created on: Aug 19, 2009
 *      Author: Adam Auton
 *      ($Revision: 230 $)
 */

#ifndef VCF_ENTRY_H_
#define VCF_ENTRY_H_

#include <algorithm>
#include <iostream>
#include <iterator>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#include <cassert>

#include "output_log.h"

using namespace std;

enum Type_enum {Integer=0, Float=1, Character=2, String=3, Flag=4};

class Field_description
{
public:
	string ID;
	int N_entries;
	Type_enum Type;
	string Description;

	Field_description() : ID(""), N_entries(0), Type(Integer), Description("") {};
	~Field_description() {};
};

class vcf_entry {
public:
	vcf_entry(const unsigned int N_indv);
	vcf_entry(const unsigned int N_indv, const string &data_line);
	~vcf_entry();

	const unsigned int N_indv;

	void parse_basic_entry(bool parse_ALT=false, bool parse_FILTER=false, bool parse_INFO=false);
	void parse_full_entry(bool parse_FORMAT=true);
	void parse_genotype_entry(unsigned int indv, bool GT=false, bool GQ=false, bool DP=false, bool FT=false);
	void parse_genotype_entries(bool GT=false, bool GQ=false, bool DP=false, bool FT=false);

	void reset(const string &vcf_data_line);

	string get_CHROM() const;
	void get_CHROM(string &out) const;
	int get_POS() const;
	string get_ID() const;
	string get_REF() const;
	string get_ALT() const;
	string get_ALT_allele(int allele_num) const;
	void get_allele(int allele_num, string &out) const;
	void get_alleles_vector(vector<string> &out) const;
	string get_FILTER() const;
	void get_FILTER_vector(vector<string> &out) const;
	double get_QUAL() const;
	string get_INFO(const set<string> &INFO_to_keep) const;
	string get_INFO_value(const string &key) const;
	string get_FORMAT() const;
	void get_indv_GENOTYPE_ids(unsigned int indv, pair<int, int> &out) const;
	void get_indv_GENOTYPE_strings(unsigned int indv, pair<string, string> &out) const;
	char get_indv_PHASE(unsigned int indv) const;
	double get_indv_GQUALITY(unsigned int indv) const;
	int get_indv_DEPTH(unsigned int indv) const;
	void get_indv_GFILTER(unsigned int indv, string &out) const;
	void get_indv_GFILTER_vector(unsigned int indv, vector<string> &out) const;
	int get_indv_ploidy(unsigned int indv) const;

	bool is_SNP() const;
	bool is_biallelic_SNP() const;
	bool is_diploid(const vector<bool> &include_indv, const vector<bool> &include_genotype) const;

	void read_indv_generic_entry(unsigned int indv, const string &FORMAT_id, string &out);
	bool FORMAT_id_exists(const string &FORMAT_id);

	unsigned int get_N_alleles() const;
	unsigned int get_N_chr(const vector<bool> &include_indv, const vector<bool> &include_genotype) const;

	void set_CHROM(const string &in);
	void set_POS(const int in);
	void set_ID(const string &in);
	void set_REF(const string &in);
	void set_ALT(const string &in);
	void set_QUAL(const double in);
	void set_FILTER(const string &FILTER_str);
	void set_FORMAT(const string &in);
	void set_INFO(const string &INFO_str);

	void add_ALT_allele(const string &in);
	void add_FILTER_entry(const string &in);
	void add_FORMAT_entry(const string &in, unsigned int pos);

	void set_indv_GENOTYPE_and_PHASE(unsigned int indv, const string &in);
	void set_indv_GENOTYPE_and_PHASE(unsigned int indv, const pair<string, string> &genotype, char phase);
	void set_indv_GENOTYPE_and_PHASE(unsigned int indv, const pair<int, int> &genotype, char phase);
	void set_indv_GENOTYPE_alleles(unsigned int indv, const pair<string, string> &in);
	void set_indv_GENOTYPE_alleles(unsigned int indv, char a1, char a2);
	void set_indv_GENOTYPE_ids(unsigned int indv, const pair<int, int> &in);
	void set_indv_PHASE(unsigned int indv, char in);
	void set_indv_GQUALITY(unsigned int indv, double in);
	void set_indv_DEPTH(unsigned int indv, int in);
	void set_indv_GFILTER(unsigned int indv, const string &in);

	void add_indv_GFILTER(unsigned int indv, const string &in);

	static void add_INFO_descriptor(const string &in);
	static void add_FILTER_descriptor(const string &in);
	static void add_FORMAT_descriptor(const string &in);

	void print(ostream &out);
	void print(ostream &out, const set<string> &INFO_to_keep, bool keep_all_INFO=false);
	void print(ostream &out, const set<string> &INFO_to_keep, bool keep_all_INFO, const vector<bool> &include_indv, const vector<bool> &include_genotype);

	void filter_genotypes_by_depth(vector<bool> &include_genotype_out, int min_depth, int max_depth);
	void filter_genotypes_by_quality(vector<bool> &include_genotype_out, double min_genotype_quality);
	void filter_genotypes_by_filter_status(vector<bool> &include_genotype_out, const set<string> &filter_flags_to_remove, bool remove_all = false);

	void get_allele_counts(vector<int> &out, unsigned int &N_non_missing_chr_out, const vector<bool> &include_indv, const vector<bool> &include_genotype) const;
	void get_genotype_counts(const vector<bool> &include_indv, const vector<bool> &include_genotype, unsigned int &out_N_hom1, unsigned int &out_N_het, unsigned int &out_N_hom2) const;

	static double SNPHWE(int obs_hets, int obs_hom1, int obs_hom2);
private:
	istringstream data_stream;

	bool basic_parsed;
	bool fully_parsed;
	bool parsed_ALT;
	bool parsed_FILTER;
	bool parsed_INFO;
	bool parsed_FORMAT;

	string CHROM;
	int POS;
	string ID;
	string REF;
	vector<string> ALT;
	double QUAL;
	vector<string> FILTER;
	bool passed_filters;
	vector<pair<string, string> > INFO;
	vector<string> FORMAT;

	vector< pair<int,int> > GENOTYPE;
	vector<int> ploidy;
	vector<char>  PHASE;
	vector<double> GQUALITY;
	vector<int>   DEPTH;
	vector< vector<string> > GFILTER;

	vector<bool> parsed_GT;
	vector<bool> parsed_GQ;
	vector<bool> parsed_DP;
	vector<bool> parsed_FT;

	int GT_idx;
	int GQ_idx;
	int DP_idx;
	int FT_idx;

	string ALT_str, FILTER_str, INFO_str, FORMAT_str, QUAL_str;
	vector<string> GENOTYPE_str;

	map<string, unsigned int> FORMAT_to_idx;

	static map<string, Field_description> INFO_map;
	static map<string, string> FILTER_map;
	static map<string, Field_description> FORMAT_map;

	static int str2int(const string &in, const int missing_value=-1);
	static double str2double(const string &in, const double missing_value=-1.0);

	static string int2str(const int in, const int missing_value=-1);
	static string double2str(const double in, const double missing_value=-1.0);

	static void tokenize(const string &in, char token, vector<string> &out);
};

#endif /* VCF_ENTRY_H_ */
