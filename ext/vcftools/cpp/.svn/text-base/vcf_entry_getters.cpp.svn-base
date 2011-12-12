/*
 * vcf_entry_getters.cpp
 *
 *  Created on: Nov 11, 2009
 *      Author: Adam Auton
 *      ($Revision: 230 $)
 */

#include "vcf_entry.h"

// Return the CHROMosome name
string vcf_entry::get_CHROM() const
{
	return CHROM;
}

// Return the CHROMosome name
void vcf_entry::get_CHROM(string &out) const
{
	out = CHROM;
}

int vcf_entry::get_POS() const
{
	return POS;
}

string vcf_entry::get_ID() const
{
	if (ID.size() == 0)
		return ".";
	return ID;
}

string vcf_entry::get_REF() const
{
	return REF;
}

string vcf_entry::get_ALT() const
{
	assert(parsed_ALT == true);

	string out;
	if (ALT.size() == 0)
		out = ".";
	else
	{
		out = ALT[0];
		for (unsigned int ui=1; ui<ALT.size(); ui++)
			out += "," + ALT[ui];
	}
	return out;
}

bool vcf_entry::is_SNP() const
{
	assert(parsed_ALT == true);

	if (REF.size() != 1)
		return false;	// Reference isn't a single base

	if (ALT.size() == 0)
		return false;	// No alternative allele

	for (unsigned int ui=0; ui<ALT.size(); ui++)
		if (ALT[ui].size() != 1)
			return false;	// Alternative allele isn't a single base

	return true;
}

bool vcf_entry::is_biallelic_SNP() const
{
	assert(parsed_ALT == true);

	if (REF.size() != 1)
		return false;	// Reference isn't a single base

	if (ALT.size() != 1)
		return false;	// Not biallelic

	if (ALT[0].size() != 1)
		return false;	// Alternative allele isn't a single base

	return true;
}

bool vcf_entry::is_diploid(const vector<bool> &include_indv, const vector<bool> &include_genotype) const
{
	for (unsigned int ui=0; ui<N_indv; ui++)
	{
		if ((include_indv[ui] == true) && (include_genotype[ui] == true))
		{
			assert(parsed_GT[ui] == true);
			if (ploidy[ui] != 2)
				return false;
		}
	}
	return true;
}

void vcf_entry::get_allele(int allele_num, string &out) const
{
	assert(parsed_ALT == true);

	if (allele_num == 0)
		out = REF;
	else if ((allele_num < 0) || (unsigned(allele_num - 1) >=  ALT.size()))
		out = ".";
	else
		out = ALT[allele_num-1];
}

string vcf_entry::get_ALT_allele(int allele_num) const
{
	assert(parsed_ALT == true);

	if ((allele_num < 0) || (unsigned(allele_num) >=  ALT.size()))
		return ".";
	return ALT[allele_num];
}

void vcf_entry::get_alleles_vector(vector<string> &out) const
{
	assert(parsed_ALT == true);
	out.resize(ALT.size()+1);
	out[0] = REF;
	copy(ALT.begin(), ALT.end(), out.begin()+1);
}

double vcf_entry::get_QUAL() const
{
	return QUAL;
}


string vcf_entry::get_FILTER() const
{
	assert(parsed_FILTER == true);

	ostringstream out;
	if ((passed_filters == false) && (FILTER.size() == 0))
		out << ".";
	else if (passed_filters == true)
		out << "PASS";
	else
	{
		out << FILTER[0];
		for (unsigned int ui=1; ui<FILTER.size(); ui++)
		{
			out << "," << FILTER[ui];
		}
	}
	return out.str();
}

void vcf_entry::get_FILTER_vector(vector<string> &out) const
{
	assert(parsed_FILTER == true);
	out = FILTER;
}


string vcf_entry::get_INFO(const set<string> &INFO_to_keep) const
{
	assert(parsed_INFO == true);

	ostringstream out;
	bool first=true;
	if ((INFO.size() > 0) && (INFO_to_keep.size() > 0))
	{
		string key;
		for (unsigned int ui=0; ui<INFO.size();ui++)
		{
			key = INFO[ui].first;
			if (INFO_to_keep.find(key) != INFO_to_keep.end())
			{
				if (first != true)
					out << ";";
				out << key << "=" << INFO[ui].second;
				first = false;
			}
		}
	}

	if (first == true)
	{	// Didn't find any INFO fields to keep
		out.str(".");
	}
	return out.str();
}

string vcf_entry::get_INFO_value(const string &key) const
{
	assert(parsed_INFO == true);

	for (unsigned int ui=0; ui<INFO.size(); ui++)
	{
		if (INFO[ui].first == key)
			return INFO[ui].second;
	}
	return "?";
}

string vcf_entry::get_FORMAT() const
{
	assert(parsed_FORMAT == true);

	string out;
	bool first = true;
	for (unsigned int ui=0; ui<FORMAT.size(); ui++)
	{
		if (first == false)
			out += ":";
		out += FORMAT[ui];
		first = false;
	}
	return out;
}

// Return the alleles of a genotype as a pair of strings.
void vcf_entry::get_indv_GENOTYPE_strings(unsigned int indv, pair<string, string> &out) const
{
	assert(parsed_GT[indv] == true);

	static string out_allele1, out_allele2;

	get_allele(GENOTYPE[indv].first, out_allele1);
	get_allele(GENOTYPE[indv].second, out_allele2);
	out = make_pair(out_allele1, out_allele2);
}


void vcf_entry::get_indv_GENOTYPE_ids(unsigned int indv, pair<int, int> &out) const
{
	assert(parsed_GT[indv] == true);
	out = GENOTYPE[indv];
}

char vcf_entry::get_indv_PHASE(unsigned int indv) const
{
	assert(parsed_GT[indv] == true);
	return PHASE[indv];
}

int vcf_entry::get_indv_DEPTH(unsigned int indv) const
{
	assert(parsed_DP[indv] == true);
	if (DEPTH.size() == 0)
		return -1;
	return DEPTH[indv];
}

double vcf_entry::get_indv_GQUALITY(unsigned int indv) const
{
	assert(parsed_GQ[indv] == true);
	if (GQUALITY.size() == 0)
		return -1;
	return GQUALITY[indv];
}

void vcf_entry::get_indv_GFILTER_vector(unsigned int indv, vector<string> &out) const
{
	assert(parsed_FT[indv] == true);
	if (GFILTER.size() > 0)
		out = GFILTER[indv];
	else
		out.resize(0);
}

void vcf_entry::get_indv_GFILTER(unsigned int indv, string &out) const
{
	assert(parsed_FT[indv] == true);

	if ((GFILTER.size() > 0) && (GFILTER[indv].size()>0))
	{
		out="";
		for (unsigned int ui=0; ui<GFILTER[indv].size(); ui++)
		{
			if (ui!=0)
				out += ";";
			out += GFILTER[indv][ui];
		}
	}
	else
		out = ".";
}

int vcf_entry::get_indv_ploidy(unsigned int indv) const
{
	assert (parsed_GT[indv]==true);
	return ploidy[indv];
}

void vcf_entry::read_indv_generic_entry(unsigned int indv, const string &FORMAT_id, string &out)
{
	if (fully_parsed == false)
		parse_full_entry(true);

	if (parsed_FORMAT == false)
		set_FORMAT(FORMAT_str);

	out = ".";

	if (FORMAT_to_idx.find(FORMAT_id) != FORMAT_to_idx.end())
	{
		unsigned int idx = FORMAT_to_idx[FORMAT_id];
		static string tmpstr;
		static istringstream ss;
		ss.clear();
		ss.str(GENOTYPE_str[indv]);

		for (unsigned int ui=0; ui <= idx; ui++)
		{
			getline(ss, tmpstr, ':');
			if (ui == idx)
			{
				out = tmpstr;
				break;
			}
			if (!ss.good())
				break;
		}
	}
}

bool vcf_entry::FORMAT_id_exists(const string &FORMAT_id)
{
	assert(parsed_FORMAT == true);
	if (FORMAT_to_idx.find(FORMAT_id) != FORMAT_to_idx.end())
		return true;
	return false;
}

unsigned int vcf_entry::get_N_alleles() const
{
	assert(parsed_ALT == true);
	return (ALT.size()+1);
}

unsigned int vcf_entry::get_N_chr(const vector<bool> &include_indv, const vector<bool> &include_genotype) const
{
	unsigned int out=0;

	for (unsigned int ui=0; ui<N_indv; ui++)
	{
		if ((include_indv[ui] == true) && (include_genotype[ui] == true))
		{
			assert(parsed_GT[ui] == true);
			out += ploidy[ui];
		}
	}
	return out;
}


// Return the frequency (counts) of each allele.
void vcf_entry::get_allele_counts(vector<int> &out, unsigned int &N_non_missing_chr_out, const vector<bool> &include_indv, const vector<bool> &include_genotype) const
{
	pair<int,int> genotype;
	vector<int> allele_counts(get_N_alleles(), 0);
	N_non_missing_chr_out = 0;

	for (unsigned int ui=0; ui<N_indv; ui++)
	{
		if ((include_indv[ui] == true) && (include_genotype[ui] == true))
		{
			assert(parsed_GT[ui] == true);
			get_indv_GENOTYPE_ids(ui, genotype);
			if (genotype.first != -1)
			{
				allele_counts[genotype.first]++;
				N_non_missing_chr_out++;
			}
			if (genotype.second != -1)
			{
				allele_counts[genotype.second]++;
				N_non_missing_chr_out++;
			}
		}
	}
	out = allele_counts;
}

// Return the counts of homozygote1, heterozygotes, and homozygote2
void vcf_entry::get_genotype_counts(const vector<bool> &include_indv, const vector<bool> &include_genotype, unsigned int &out_N_hom1, unsigned int &out_N_het, unsigned int &out_N_hom2) const
{
	out_N_hom1 = 0; out_N_hom2 = 0; out_N_het = 0;
	pair<int, int> genotype;
	if (ALT.size() > 1)
		error("Tried to return the genotype counts of a non-biallelic SNP", 99);

	for (unsigned int ui=0; ui<N_indv; ui++)
	{
		if ((include_indv[ui] == true) && (include_genotype[ui] == true))
		{
			assert(parsed_GT[ui] == true);
			get_indv_GENOTYPE_ids(ui, genotype);
			if ((genotype.first != -1) && (genotype.second != -1))
			{
				if (genotype.first != genotype.second)
				{
					out_N_het++;
				}
				else if (genotype.first == 0)
				{
					out_N_hom1++;
				}
				else if (genotype.first == 1)
				{
					out_N_hom2++;
				}
				else
				{
					error("Unknown allele in genotype", 98);
				}
			}
		}
	}
}
