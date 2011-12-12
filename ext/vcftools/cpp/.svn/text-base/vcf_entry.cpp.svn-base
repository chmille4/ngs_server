/*
 * vcf_entry.cpp
 *
 *  Created on: Aug 19, 2009
 *      Author: Adam Auton
 *      ($Revision: 230 $)
 */

#include "vcf_entry.h"

map<string, Field_description> vcf_entry::INFO_map;
map<string, string> vcf_entry::FILTER_map;
map<string, Field_description> vcf_entry::FORMAT_map;

// Create a VCF on the basis of a data line.
vcf_entry::vcf_entry(const unsigned int N_indv, const string &line)
	:  N_indv(N_indv),
	   data_stream(line),
	   basic_parsed(false), fully_parsed(false),
	   parsed_ALT(false), parsed_FILTER(false),
	   parsed_INFO(false), parsed_FORMAT(false),
	   CHROM(""), POS(-1), REF(""), QUAL(-1),
	   passed_filters(false),
	   parsed_GT(N_indv, false), parsed_GQ(N_indv, false), parsed_DP(N_indv, false),
	   parsed_FT(N_indv, false),
	   ALT_str(""), FILTER_str(""), INFO_str(""), FORMAT_str(""), QUAL_str("")
{
}

// Create an empty VCF entry
vcf_entry::vcf_entry(const unsigned int N_indv)
	: N_indv(N_indv),
	  data_stream("0\t0\t.\tN\t.\t.\t.\t."),
	  basic_parsed(false), fully_parsed(false),
	  parsed_ALT(false), parsed_FILTER(false),
	  parsed_INFO(false), parsed_FORMAT(false),
	  CHROM(""), POS(-1), REF(""), QUAL(-1),
	  passed_filters(false),
	  parsed_GT(N_indv, false), parsed_GQ(N_indv, false), parsed_DP(N_indv, false),
	  parsed_FT(N_indv, false),
	  ALT_str(""), FILTER_str(""), INFO_str(""), FORMAT_str(""), QUAL_str("")
{
}

vcf_entry::~vcf_entry() {}

// Reset the VCF entry object with a new data line
void vcf_entry::reset(const string &vcf_data_line)
{
	basic_parsed = false;
	fully_parsed = false;
	parsed_ALT = false;
	parsed_FILTER = false;
	parsed_INFO = false;
	parsed_FORMAT = false;

	data_stream.clear();
	data_stream.str(vcf_data_line);

	fill(parsed_GT.begin(), parsed_GT.end(), false);
	fill(parsed_GQ.begin(), parsed_GQ.end(), false);
	fill(parsed_DP.begin(), parsed_DP.end(), false);
	fill(parsed_FT.begin(), parsed_FT.end(), false);
}

// Tokenize the basic information in a VCF data line (at the tab level)
void vcf_entry::parse_basic_entry(bool parse_ALT, bool parse_FILTER, bool parse_INFO)
{
	data_stream >> CHROM >> POS >> ID >> REF >> ALT_str >> QUAL_str >> FILTER_str >> INFO_str;
	QUAL = str2double(QUAL_str);

	// Convert to uppercase for consistency
	// Note that VCF v4.1 allows mixtures of lower/upper case in REF and ALT.
	// However, the spec specifically states that tools using VCF are not required
	// to preserve the case.
	std::transform(REF.begin(), REF.end(), REF.begin(), ::toupper);
	std::transform(ALT_str.begin(), ALT_str.end(),ALT_str.begin(), ::toupper);

	parsed_ALT = false;
	parsed_FILTER = false;
	parsed_INFO = false;
	basic_parsed = true;

	if (parse_ALT)
		set_ALT(ALT_str);
	if (parse_FILTER)
		set_FILTER(FILTER_str);
	if (parse_INFO)
		set_INFO(INFO_str);
}

// Tokenize the genotype information (at the 'tab' level) in the VCF entry
void vcf_entry::parse_full_entry(bool parse_FORMAT)
{
	if (basic_parsed == false)
		parse_basic_entry();

	data_stream >> FORMAT_str;

	if (parse_FORMAT)
		set_FORMAT(FORMAT_str);

	string tmpstr; tmpstr.reserve(64);
	GENOTYPE_str.resize(N_indv, tmpstr);

	for (unsigned int ui=0; ui<N_indv; ui++)
		data_stream >> GENOTYPE_str[ui];

	// The following line copies the GENOTYPE fields from the stringstream into the GENOTYPE_str vector.
	// Is actually slower than the above code.
	//copy(istream_iterator<string>(data_stream), istream_iterator<string>(), GENOTYPE_str.begin());

	fully_parsed = true;
}

// Tokenize a given genotype entry into it's component parts
void vcf_entry::parse_genotype_entry(unsigned int indv, bool GT, bool GQ, bool DP, bool FT)
{
	if (fully_parsed == false)
		parse_full_entry(true);

	if (parsed_FORMAT == false)
		set_FORMAT(FORMAT_str);

	static string tmpstr;
	static istringstream ss;
	ss.clear(); ss.str(GENOTYPE_str[indv]);

	int N_required = GT + GQ + DP + FT;
	int N_got = 0;

	int i=0;
	while (getline(ss, tmpstr, ':'))
	{
		if (GT && (i == GT_idx)) // (FORMAT[ui] == "GT")
		{
			set_indv_GENOTYPE_and_PHASE(indv, tmpstr);
			N_got++;
		}
		else if (GQ && (i == GQ_idx)) // (FORMAT[ui] == "GQ")
		{
			set_indv_GQUALITY(indv, str2double(tmpstr));
			N_got++;
		}
		else if (DP && (i == DP_idx)) // (FORMAT[ui] == "DP")
		{
			set_indv_DEPTH(indv, str2int(tmpstr));
			N_got++;
		}
		else if (FT && (i == FT_idx)) // (FORMAT[ui] == "FT")
		{
			set_indv_GFILTER(indv, tmpstr);
			N_got++;
		}

		if (N_got == N_required)
			break;
		i++;
	}

	// Set missing return values if requested a value, but couldn't find it
	if (GT && (parsed_GT[indv] == false))
	{
		set_indv_GENOTYPE_and_PHASE(indv, make_pair(-1,-1), '/');
	}
	if (GQ && (parsed_GQ[indv] == false))
	{
		set_indv_GQUALITY(indv, -1);
	}
	if (DP && (parsed_DP[indv] == false))
	{
		set_indv_DEPTH(indv, -1);
	}
	if (FT && (parsed_FT[indv] == false))
	{
		set_indv_GFILTER(indv, "");
	}
}

// Read the VCF entry and fully populate the object
void vcf_entry::parse_genotype_entries(bool GT, bool GQ, bool DP, bool FT)
{
	for (unsigned int ui=0; ui<N_indv; ui++)
		parse_genotype_entry(ui, GT, GQ, DP, FT);
}

void vcf_entry::print(ostream &out)
{
	vector<bool> include_indv(N_indv, true);
	vector<bool> include_genotype(N_indv, true);
	set<string> INFO_to_keep;
	print(out, INFO_to_keep, false, include_indv, include_genotype);
}

void vcf_entry::print(ostream &out, const set<string> &INFO_to_keep, bool keep_all_INFO)
{
	vector<bool> include_indv(N_indv, true);
	vector<bool> include_genotype(N_indv, true);
	print(out, INFO_to_keep, keep_all_INFO, include_indv, include_genotype);
}

// Output VCF entry to output stream
void vcf_entry::print(ostream &out, const set<string> &INFO_to_keep, bool keep_all_INFO, const vector<bool> &include_indv, const vector<bool> &include_genotype)
{
	if (fully_parsed == false)
		parse_full_entry();

	out << get_CHROM() << '\t' << POS << '\t' << get_ID() << '\t' << REF << '\t' << get_ALT();

	out << '\t' << double2str(QUAL);
	out << '\t' << get_FILTER();
	if (keep_all_INFO == false)
		out << '\t' << get_INFO(INFO_to_keep);
	else
		out << '\t' << INFO_str;

	pair<int, int> genotype;
	string GFILTER_tmp;
	if (FORMAT.size() > 0)
	{
		char PHASE;
		out << '\t' << get_FORMAT();

		for (unsigned int ui=0; ui<N_indv; ui++)
		{
			if (include_indv[ui] == false)
				continue;
			out << '\t';
			for (int count=0; count<(int)FORMAT.size(); count++)
			{
				if (count == GT_idx) // (FORMAT[count] == "GT")
				{
					if (count != 0)	out << ':';
					if (include_genotype[ui] == true)
					{
						get_indv_GENOTYPE_ids(ui, genotype);
						PHASE = get_indv_PHASE(ui);
						if ((genotype.first != -1) && (genotype.second != -1))
							out << int2str(genotype.first) << PHASE << int2str(genotype.second);
						else if ((PHASE == '|') && (genotype.second == -1))
							out << int2str(genotype.first);	// Handle haploid case
						else
							out << int2str(genotype.first) << PHASE << int2str(genotype.second);
					}
					else
						out << "./.";
				}
				else if (count == GQ_idx) //(FORMAT[count] == "GQ")
				{
					if (count != 0)	out << ':';
					out << double2str(get_indv_GQUALITY(ui));
				}
				else if (count == DP_idx) // (FORMAT[count] == "DP")
				{
					if (count != 0)	out << ':';
					out << int2str(get_indv_DEPTH(ui));
				}
				else if (count == FT_idx) // (FORMAT[count] == "FT")
				{
					if (count != 0)	out << ':';
					get_indv_GFILTER(ui, GFILTER_tmp);
					out << GFILTER_tmp;
				}
				else
				{	// Unknown FORMAT so just replicate original output
					if (count != 0)	out << ':';
					read_indv_generic_entry(ui, FORMAT[count], GFILTER_tmp);
					out << GFILTER_tmp;
				}
			}
		}
	}
//	out << endl;
	out << '\n';	// endl flushes the buffer, which is slow. This (should be) quicker.
}

// Set the include_genotype flag on the basis of depth
void vcf_entry::filter_genotypes_by_depth(vector<bool> &include_genotype_out, int min_depth, int max_depth)
{
	if (fully_parsed == false)
		parse_full_entry();

	//if (FORMAT_to_idx.find("DP") != FORMAT_to_idx.end())
	if (DP_idx != -1)
	{	// Have depth info
		int depth;
		include_genotype_out.resize(N_indv, true);
		for (unsigned int ui=0; ui<N_indv; ui++)
		{
			if (parsed_DP[ui] == false)
				parse_genotype_entry(ui, false, false, true);
			depth = get_indv_DEPTH(ui);
			if ((depth < min_depth) || (depth > max_depth))
				include_genotype_out[ui] = false;
		}
	}
}

// Filter specific genotypes by quality
void vcf_entry::filter_genotypes_by_quality(vector<bool> &include_genotype_out, double min_genotype_quality)
{
	if (fully_parsed == false)
		parse_full_entry();

	//if (FORMAT_to_idx.find("GQ") != FORMAT_to_idx.end())
	if (GQ_idx != -1)
	{	// Have quality info
		double quality;
		include_genotype_out.resize(N_indv, true);
		for (unsigned int ui=0; ui<N_indv; ui++)
		{
			if (parsed_GQ[ui] == false)
				parse_genotype_entry(ui, false, true);
			quality = get_indv_GQUALITY(ui);
			if (quality < min_genotype_quality)
				include_genotype_out[ui] = false;
		}
	}
}

// Exclude genotypes with a filter flag.
void vcf_entry::filter_genotypes_by_filter_status(vector<bool> &include_genotype_out, const set<string> &filter_flags_to_remove, bool remove_all)
{
	if (fully_parsed == false)
		parse_full_entry();

	string filter;
	vector<string> GFILTERs;
	//if (FORMAT_to_idx.find("FT") != FORMAT_to_idx.end())
	if (FT_idx != -1)
	{	// Have GFilter info
		include_genotype_out.resize(N_indv, true);
		for (unsigned int ui=0; ui<N_indv; ui++)
		{
			if (parsed_FT[ui] == false)
				parse_genotype_entry(ui, false, false, false, true);
			get_indv_GFILTER_vector(ui, GFILTERs);

			if ((remove_all == true) && (GFILTERs.size() > 0))
				include_genotype_out[ui] = false;
			else
			{
				for (unsigned int uj=0; uj<GFILTERs.size(); uj++)
					if (filter_flags_to_remove.find(GFILTERs[uj]) != filter_flags_to_remove.end())
						include_genotype_out[ui] = false;
			}
		}
	}
}

/*
// This function implements an exact SNP test of Hardy-Weinberg
// Equilibrium as described in Wigginton, JE, Cutler, DJ, and
// Abecasis, GR (2005) A Note on Exact Tests of Hardy-Weinberg
// Equilibrium. American Journal of Human Genetics. 76: 000 - 000
//
// Written by Jan Wigginton
*/
double vcf_entry::SNPHWE(int obs_hets, int obs_hom1, int obs_hom2)
{
	if (obs_hom1 + obs_hom2 + obs_hets == 0 ) return 1;

	if (obs_hom1 < 0 || obs_hom2 < 0 || obs_hets < 0)
		error("Internal error: negative count in HWE test", 91);

	int obs_homc = obs_hom1 < obs_hom2 ? obs_hom2 : obs_hom1;
	int obs_homr = obs_hom1 < obs_hom2 ? obs_hom1 : obs_hom2;

	int rare_copies = 2 * obs_homr + obs_hets;
	int genotypes   = obs_hets + obs_homc + obs_homr;

	double * het_probs = (double *) malloc((size_t) (rare_copies + 1) * sizeof(double));
	if (het_probs == NULL)
		error("Internal error: SNP-HWE: Unable to allocate array", 90);

	for (int i = 0; i <= rare_copies; i++)
		het_probs[i] = 0.0;

	/* start at midpoint */
	int mid = rare_copies * (2 * genotypes - rare_copies) / (2 * genotypes);

	/* check to ensure that midpoint and rare alleles have same parity */
	if ((rare_copies & 1) ^ (mid & 1))
		mid++;

	int curr_hets = mid;
	int curr_homr = (rare_copies - mid) / 2;
	int curr_homc = genotypes - curr_hets - curr_homr;

	het_probs[mid] = 1.0;
	double sum = het_probs[mid];
	for (curr_hets = mid; curr_hets > 1; curr_hets -= 2)
	{
	  het_probs[curr_hets - 2] = het_probs[curr_hets] * curr_hets * (curr_hets - 1.0) / (4.0 * (curr_homr + 1.0) * (curr_homc + 1.0));
	  sum += het_probs[curr_hets - 2];

	  /* 2 fewer heterozygotes for next iteration -> add one rare, one common homozygote */
	  curr_homr++;
	  curr_homc++;
	}

	curr_hets = mid;
	curr_homr = (rare_copies - mid) / 2;
	curr_homc = genotypes - curr_hets - curr_homr;
	for (curr_hets = mid; curr_hets <= rare_copies - 2; curr_hets += 2)
	{
	  het_probs[curr_hets + 2] = het_probs[curr_hets] * 4.0 * curr_homr * curr_homc	/((curr_hets + 2.0) * (curr_hets + 1.0));
	  sum += het_probs[curr_hets + 2];

	  /* add 2 heterozygotes for next iteration -> subtract one rare, one common homozygote */
	  curr_homr--;
	  curr_homc--;
	}

	for (int i = 0; i <= rare_copies; i++)
		het_probs[i] /= sum;

	/* alternate p-value calculation for p_hi/p_lo
	double p_hi = het_probs[obs_hets];
	for (int i = obs_hets + 1; i <= rare_copies; i++)
	 p_hi += het_probs[i];

	double p_lo = het_probs[obs_hets];
	for (int i = obs_hets - 1; i >= 0; i--)
	  p_lo += het_probs[i];

	double p_hi_lo = p_hi < p_lo ? 2.0 * p_hi : 2.0 * p_lo;
	*/

	double p_hwe = 0.0;
	/*  p-value calculation for p_hwe  */
	for (int i = 0; i <= rare_copies; i++)
	{
		if (het_probs[i] > het_probs[obs_hets])
			continue;
		p_hwe += het_probs[i];
	}

	p_hwe = p_hwe > 1.0 ? 1.0 : p_hwe;

	free(het_probs);

	return p_hwe;
}

int vcf_entry::str2int(const string &in, const int missing_value)
{
	if ((in.size() == 0) || (in == "."))
		return missing_value;
	else
		return atoi(in.c_str());
}

double vcf_entry::str2double(const string &in, const double missing_value)
{
	if ((in.size() == 0) || (in == "."))
		return missing_value;
	else
		return atof(in.c_str());
}

string vcf_entry::int2str(const int in, const int missing_value)
{
	if (in == missing_value)
		return ".";
	else
	{
		static ostringstream out;
		out.str(""); out.clear();
		out << in;
		return out.str();
	}
}

string vcf_entry::double2str(const double in, const double missing_value)
{
	if (in == missing_value)
		return ".";
	else
	{
		static ostringstream out;
		out.str(""); out.clear();
		out << in;
		return out.str();
	}
}

void vcf_entry::tokenize(const string &in, char token, vector<string> &out)
{
	out.resize(0);
	istringstream ss(in);
	string tmp;
	while( getline(ss, tmp, token) )
	{
		out.push_back(tmp);
	}
}
