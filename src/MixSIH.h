#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <vector>
#include <string>
#include <algorithm>
#include <cassert>

using namespace std;

double digamma(double x){
	if(x <=  0){
		printf("error at digamma function :%lf\n", x);
	}
	
	return log(x) - 1.0/(2*x) - 1.0/(12*x*x) + 1.0/(120*x*x*x*x) - 1/(252*x*x*x*x*x*x);
}

//Fragment class stores information of one fragment.
class Fragment{
private:
	int _size;
	vector<int> _pos;
	vector<char> _base;

public:
	Fragment(){
		_size = 0;
	}
	
	void Add(int pos, char base){
		_pos.push_back(pos);
		_base.push_back(base);
		_size++;
	}
	
	void Offset(int offset){
		for(int i=0; i<_size; i++){
			_pos[i] -= offset;
		}
	}
	
	int Size(){
		return _size;
	}
	
	int Pos(int num){
		return _pos[num];
	}
	
	char Base(int num){
		return _base[num];
	}
};

class SNP{
private:
	int _size;
	int _cont;
	double _connectivity;
	//fragments num which contain this SNP
	vector<int> _fragments;

public:	
	SNP(){
		_size = 0;
		_cont = 0;
		_connectivity = 0;
	}
	
	void Assign_fragment(int num){
		_fragments.push_back(num);
		_size++;
	}
	
	int Size(){
		return _size;
	}
	
	int Fragments(int num){
		return _fragments[num];
	}
	
	void Connectivity(double num){
		_connectivity = num;
	}
	
	double Connectivity(){
		return _connectivity;
	}
	
	void Add_cont(){
		_cont++;
	}
	
	int Cont(){
		return _cont;
	}
};

class Haplotype{
private:
	int _snp_num, _frag_num, _offset;
	double _alpha, _lambda0;
	vector<vector<double> > _lambda;
	vector<vector<double> > _old_lambda;
	vector<vector<double> > _gamma1;
	vector<vector<vector<vector<double> > > > _gamma2;
	//_fragments stores the 
	vector<Fragment> _fragments;
	vector<SNP> _snps;
	vector<int> _flag;

public:	
	void Set_fragment(vector<Fragment> &all_fragments, vector<int> &id, vector<int> &pos, int offset){
		_fragments.resize(id.size());
		for(unsigned int i=0; i<id.size(); i++){
			_fragments[i] = all_fragments[id[i]];
		}
		
		_frag_num = _fragments.size();
		_snp_num = pos[pos.size()-1] - pos[0] + 1;
		_offset = offset;
	}
	
	void Set_fragment(vector<Fragment> &fragments, int num_of_SNP, int offset){
		_fragments.resize(fragments.size());
		for(unsigned int i=0; i<fragments.size(); i++){
			_fragments[i] = fragments[i];
		}
		
		_frag_num = fragments.size();
		_snp_num = num_of_SNP;
		_offset = offset;
	}
	
	//SNPs[i] stores id of the fragments which cover the i-th SNP
	void Assign_fragment_to_SNP(){
		int start, end;
		_snps.resize(_snp_num);
		
		for(int i=0; i<_frag_num; i++){
			//j start from 1
			start = _fragments[i].Pos(0)+1;
			end = _fragments[i].Pos(_fragments[i].Size()-1);
			
			for(int j=start; j<=end; j++){
				_snps[j].Assign_fragment(i);
			}
		}
	}

	void Init(double alpha){
		_lambda0 = 0.5;
		_alpha = alpha;
		
		srand(time(NULL));
		_lambda.resize(_snp_num, vector<double>(2, _lambda0));
		_old_lambda.resize(_snp_num, vector<double>(2, _lambda0));
		//initialize lambda
		for(int j=0; j<_snp_num; j++){
			for(int k=0; k<2; k++){
				_lambda[j][k] += 0.1*(double)rand()/RAND_MAX;
			}
		}
		
		_gamma1.resize(_frag_num, vector<double>(2, 0));
		_gamma2.resize(_frag_num, vector<vector<vector<double> > >(2, vector<vector<double> >()));
		for(int i=0; i<_frag_num; i++){
			for(int j=0; j<2; j++){
				_gamma2[i][j].resize(_fragments[i].Size(), vector<double>(2, 0));
			}
		}
		
		_flag.resize(_snp_num, 0);
	}
	
	//The probability of a fragment with normal parameter.
	double Probability1(int frag_id, int hap_id){
		double ret = 1;
		if(hap_id == 0){
			for(int i=0; i<_fragments[frag_id].Size(); i++){
				if(_fragments[frag_id].Base(i) == 0){
					ret *= ((1-_alpha)*_lambda[_fragments[frag_id].Pos(i)][0] + _alpha*_lambda[_fragments[frag_id].Pos(i)][1])/(_lambda[_fragments[frag_id].Pos(i)][0] + _lambda[_fragments[frag_id].Pos(i)][1]);
				}
				else{
					ret *= ((1-_alpha)*_lambda[_fragments[frag_id].Pos(i)][1] + _alpha*_lambda[_fragments[frag_id].Pos(i)][0])/(_lambda[_fragments[frag_id].Pos(i)][0] + _lambda[_fragments[frag_id].Pos(i)][1]);
				}
			}
		}
		else{
			for(int i=0; i<_fragments[frag_id].Size(); i++){
				if(_fragments[frag_id].Base(i) == 0){
					ret *= ((1-_alpha)*_lambda[_fragments[frag_id].Pos(i)][1] + _alpha*_lambda[_fragments[frag_id].Pos(i)][0])/(_lambda[_fragments[frag_id].Pos(i)][0] + _lambda[_fragments[frag_id].Pos(i)][1]);
				}
				else{
					ret *= ((1-_alpha)*_lambda[_fragments[frag_id].Pos(i)][0] + _alpha*_lambda[_fragments[frag_id].Pos(i)][1])/(_lambda[_fragments[frag_id].Pos(i)][0] + _lambda[_fragments[frag_id].Pos(i)][1]);
				}
			}
		}

		return ret;
	}

	//The probability of a fragment with switched parameter.
	double Probability2(int frag_id, int hap_id, int switch_pos){
		double ret = 1;
		if(hap_id == 0){
			for(int i=0; i<_fragments[frag_id].Size(); i++){
				if(_fragments[frag_id].Pos(i) <= switch_pos){
					if(_fragments[frag_id].Base(i) == 0){
						ret *= ((1-_alpha)*_lambda[_fragments[frag_id].Pos(i)][0] + _alpha*_lambda[_fragments[frag_id].Pos(i)][1])/(_lambda[_fragments[frag_id].Pos(i)][0] + _lambda[_fragments[frag_id].Pos(i)][1]);
					}
					else{
						ret *= ((1-_alpha)*_lambda[_fragments[frag_id].Pos(i)][1] + _alpha*_lambda[_fragments[frag_id].Pos(i)][0])/(_lambda[_fragments[frag_id].Pos(i)][0] + _lambda[_fragments[frag_id].Pos(i)][1]);
					}
				}
				else{
					if(_fragments[frag_id].Base(i) == 0){
						ret *= ((1-_alpha)*_lambda[_fragments[frag_id].Pos(i)][1] + _alpha*_lambda[_fragments[frag_id].Pos(i)][0])/(_lambda[_fragments[frag_id].Pos(i)][0] + _lambda[_fragments[frag_id].Pos(i)][1]);
					}
					else{
						ret *= ((1-_alpha)*_lambda[_fragments[frag_id].Pos(i)][0] + _alpha*_lambda[_fragments[frag_id].Pos(i)][1])/(_lambda[_fragments[frag_id].Pos(i)][0] + _lambda[_fragments[frag_id].Pos(i)][1]);
					}
				}
			}
		}
		else{
			for(int i=0; i<_fragments[frag_id].Size(); i++){
				if(_fragments[frag_id].Pos(i) <= switch_pos){
					if(_fragments[frag_id].Base(i) == 0){
						ret *= ((1-_alpha)*_lambda[_fragments[frag_id].Pos(i)][1] + _alpha*_lambda[_fragments[frag_id].Pos(i)][0])/(_lambda[_fragments[frag_id].Pos(i)][0] + _lambda[_fragments[frag_id].Pos(i)][1]);
					}
					else{
						ret *= ((1-_alpha)*_lambda[_fragments[frag_id].Pos(i)][0] + _alpha*_lambda[_fragments[frag_id].Pos(i)][1])/(_lambda[_fragments[frag_id].Pos(i)][0] + _lambda[_fragments[frag_id].Pos(i)][1]);
					}
				}
				else{
					if(_fragments[frag_id].Base(i) == 0){
						ret *= ((1-_alpha)*_lambda[_fragments[frag_id].Pos(i)][0] + _alpha*_lambda[_fragments[frag_id].Pos(i)][1])/(_lambda[_fragments[frag_id].Pos(i)][0] + _lambda[_fragments[frag_id].Pos(i)][1]);
					}
					else{
						ret *= ((1-_alpha)*_lambda[_fragments[frag_id].Pos(i)][1] + _alpha*_lambda[_fragments[frag_id].Pos(i)][0])/(_lambda[_fragments[frag_id].Pos(i)][0] + _lambda[_fragments[frag_id].Pos(i)][1]);
					}
				}
			}
		}
		return ret;
	}
	
	double Log_likelihood(){
		double ret=0, tmp;
		
		for(int i=0; i<_frag_num; i++){
			tmp = 0;
			for(int j=0; j<2; j++){
				tmp += 0.5 * Probability1(i, j);
			}
			ret += log(tmp);
		}
		
		return ret;
	}

	void E_step(){
		double sum=0, digamma_lambda;

		for(int i=0; i<_frag_num; i++){		
			for(int j=0; j<_fragments[i].Size(); j++){
				digamma_lambda = digamma(_lambda[_fragments[i].Pos(j)][0] + _lambda[_fragments[i].Pos(j)][1]);

				if(_fragments[i].Base(j) == 0){
					_gamma2[i][0][j][0] = exp(log(1-_alpha) + digamma(_lambda[_fragments[i].Pos(j)][0]) - digamma_lambda);
					_gamma2[i][0][j][1] = exp(log(_alpha) + digamma(_lambda[_fragments[i].Pos(j)][1]) - digamma_lambda);
					_gamma2[i][1][j][0] = exp(log(1-_alpha) + digamma(_lambda[_fragments[i].Pos(j)][1]) - digamma_lambda);
					_gamma2[i][1][j][1] = exp(log(_alpha) + digamma(_lambda[_fragments[i].Pos(j)][0]) - digamma_lambda);
				}
				else{
					_gamma2[i][0][j][0] = exp(log(1-_alpha) + digamma(_lambda[_fragments[i].Pos(j)][1]) - digamma_lambda);
					_gamma2[i][0][j][1] = exp(log(_alpha) + digamma(_lambda[_fragments[i].Pos(j)][0]) - digamma_lambda);
					_gamma2[i][1][j][0] = exp(log(1-_alpha) + digamma(_lambda[_fragments[i].Pos(j)][0]) - digamma_lambda);
					_gamma2[i][1][j][1] = exp(log(_alpha) + digamma(_lambda[_fragments[i].Pos(j)][1]) - digamma_lambda);
				}
			}
			
			double tmp = 1.0;
			for(int j=0; j<_fragments[i].Size(); j++){
				tmp *= (_gamma2[i][1][j][0] + _gamma2[i][1][j][1])/(_gamma2[i][0][j][0] + _gamma2[i][0][j][1]);
			}
			for(int j=0; j<2; j++){
				if(isinf(tmp)){
					if(j == 0){
						_gamma1[i][j] = 0.0;
					}
					else{
						_gamma1[i][j] = 1.0;
					}
				}				
				else{
					if(j == 0){
						_gamma1[i][j] = 1.0 / (1.0 + tmp);
					}
					else{
						_gamma1[i][j] = tmp / (1.0 + tmp);
					}
				}
			}
		}

		for(int i=0; i<_frag_num; i++){
			for(int j=0; j<_fragments[i].Size(); j++){
				for(int k=0; k<2; k++){
					sum = _gamma2[i][k][j][0] + _gamma2[i][k][j][1];
					_gamma2[i][k][j][0] /= sum;
					_gamma2[i][k][j][1] /= sum;
				}
			}
		}	
	}
	
	void M_step(){
		//update lambda
		for(int i=0; i<_snp_num; i++){
			for(int j=0; j<2; j++){
				_old_lambda[i][j] = _lambda[i][j];
				_lambda[i][j] = _lambda0;
			}
		}
		for(int i=0; i<_frag_num; i++){
			for(int k=0; k<_fragments[i].Size(); k++){
				if(_fragments[i].Base(k) == 0){
					_lambda[_fragments[i].Pos(k)][0] += _gamma1[i][0]*_gamma2[i][0][k][0] + _gamma1[i][1]*_gamma2[i][1][k][1];
					_lambda[_fragments[i].Pos(k)][1] += _gamma1[i][0]*_gamma2[i][0][k][1] + _gamma1[i][1]*_gamma2[i][1][k][0];
				}
				else{
					_lambda[_fragments[i].Pos(k)][0] += _gamma1[i][0]*_gamma2[i][0][k][1] + _gamma1[i][1]*_gamma2[i][1][k][0];
					_lambda[_fragments[i].Pos(k)][1] += _gamma1[i][0]*_gamma2[i][0][k][0] + _gamma1[i][1]*_gamma2[i][1][k][1];
				}
			}
		}
	}
	
	bool check_converge(){
		for(int i=0; i<_snp_num; i++){
			for(int j=0; j<2; j++){
				if(fabs(_lambda[i][j] - _old_lambda[i][j]) > 0.0001)
					return 1;
			}
		}
		
		return 0;
	}
	
	void EM(){
		int change_iteration=_snp_num;
		int em_iteration=100, min_pos, flag=0;
		double min, ll, changed_ll;
		Haplotype changed_hap;
		changed_hap.Set_fragment(_fragments, _snp_num, _offset);
		changed_hap.Assign_fragment_to_SNP();
		changed_hap.Init(_alpha);
		
		for(int i=0; i<em_iteration; i++){
			E_step();
			M_step();

			if(check_converge() == 0){
				break;
			}
		}
		ll = Log_likelihood();
		Check_connectivity(0);

		for(int i=0; i<change_iteration; i++){
			if(flag==1){
				Check_connectivity(0);
			}
			min_pos = -1;
			min = 1000000;
			for(int j=1; j<_snp_num; j++){
				if(_snps[j].Connectivity()!=-1 && min > _snps[j].Connectivity() && _flag[j] < 2){
					min_pos = j;
					min = _snps[j].Connectivity();
				}
			}
			if(min_pos == -1){
				break;
			}
		
			if(min > 7.0){
				break;
			}
			
			_flag[min_pos]++;
			
			for(int j=0; j<min_pos; j++){
				for(int k=0; k<2; k++){
					changed_hap._lambda[j][k] = _lambda[j][k];
				}
			}
			
			for(int j=min_pos; j<_snp_num; j++){
				changed_hap._lambda[j][0] = _lambda[j][1];
				changed_hap._lambda[j][1] = _lambda[j][0];
			}
			
			changed_hap.EM_nochange();
			
			changed_ll = changed_hap.Log_likelihood();
			
			if(changed_ll > ll){
				flag = 1;
				ll = changed_ll;
				for(int j=0; j<_snp_num; j++){
					for(int k=0; k<2; k++)
						_lambda[j][k] = changed_hap._lambda[j][k];
				}
			}
		}
		
		Check_connectivity(1);
	}
			
	void EM_nochange(){
		int iteration=100;

		for(int i=0; i<iteration; i++){
			E_step();
			M_step();
			
			if(check_converge() == 0){
				break;
			}
		}
	}
	
	//check cennection of haplotype
	void Check_connectivity(int flag){
		for(int i=1; i<_snp_num; i++){
			if(flag == 0 && _flag[i] == 2){
				continue;
			}

			if(_snps[i].Size() == 0){
				_snps[i].Connectivity(-1);
				continue;
			}
			
			double log_p1=0, log_p2=0;
			double tmp;
			
			//The log probability with normal parameter.
			for(int j=0; j<_snps[i].Size(); j++){
				tmp = 0;
				for(int k=0; k<2; k++){
					tmp += 0.5 * Probability1(_snps[i].Fragments(j), k);
				}
				log_p1 += log(tmp);
			}
			
			//The log probability with switched parameter.
			for(int j=0; j<_snps[i].Size(); j++){
				tmp = 0.5 * Probability2(_snps[i].Fragments(j), 0, i-1);
				tmp += 0.5 * Probability2(_snps[i].Fragments(j), 1, i-1);
			
				log_p2 += log(tmp);
			}
			
			_snps[i].Connectivity(log_p1-log_p2);
		}
	}

	void Output_profile(FILE *fp){
		int cont = _snp_num;
		for(int i=0; i<_snp_num; i++){
			if(_lambda[i][0] == _lambda[i][1]){
				cont--;
			}
		}

		fprintf(fp, "BLOCK: offset: %d len: %d phased: %d\n", _offset, _snp_num, cont);

		double sum;
		for(int i=0; i<_snp_num; i++){
			fprintf(fp, "%d\t", i+_offset);
			sum = _lambda[i][0] + _lambda[i][1];
				
			fprintf(fp, "%.3lf\t%.3lf\t%.3lf\n", _lambda[i][0]/sum, _lambda[i][1]/sum, _snps[i].Connectivity());
		}
		fprintf(fp, "********\n");
	}
	
	void Output_haplotype(FILE *fp){
		int cont = _snp_num;
		for(int i=0; i<_snp_num; i++){
			if(_lambda[i][0] == _lambda[i][1]){
				cont--;
			}
		}

		fprintf(fp, "BLOCK: offset: %d len: %d phased: %d\n", _offset, _snp_num, cont);

		for(int i=0; i<_snp_num; i++){
			if(_lambda[i][0] > _lambda[i][1]){
				fprintf(fp, "%d\t0\t1\n", i+_offset);
			}
			else if(_lambda[i][0] < _lambda[i][1]){
				fprintf(fp, "%d\t1\t0\n", i+_offset);
			}
			else{
				fprintf(fp, "%d\t-\t-\n", i+_offset);
			}
		}
		fprintf(fp, "********\n");
	}
};

//position is 1-origin
//input file must be sorted by start site of fragment
void set_fragment(FILE *fp, vector<Haplotype> &hap){
	char buf[1000000];
	string tmp="";
	
	int i, all_fragment;
	fgets(buf, 1000000, fp);
	for(i=0; buf[i]!=' ' && buf[i]!='\r' && buf[i]!='\n'; i++)
		tmp += buf[i];
	all_fragment = atoi(tmp.c_str()) - 1;
	
	vector<Fragment> all_fragments(all_fragment);
	
	int pos, contig_num;
	for(int j=0; j<all_fragment; j++){
		fgets(buf, 1000000, fp);
		
		tmp = "";
		for(i=0; ; i++){
			if(buf[i] == ' '){
				contig_num = atoi(tmp.c_str());
				break;
			}
			tmp += buf[i];
		}
		
		for(i++;;i++){
			if(buf[i] == ' ')
				break;
		}
		
		for(int k=0; k<contig_num; k++){
			tmp = "";
			for(i++; buf[i]!=' '; i++){
				tmp += buf[i];
			}
			pos = atoi(tmp.c_str());

			for(i++; buf[i]!=' ' && buf[i]!='\n' && buf[i] != '\r'; i++, pos++){
				//relative position of each block
				if(buf[i] == '0'){
					all_fragments[j].Add(pos, 0);
				}
				else if(buf[i] == '1'){
					all_fragments[j].Add(pos, 1);
				}
			}
		}
	}
	
	int flag, flag2=0, max;
	vector<vector<int> > block(all_fragment);
	vector<vector<int> > block_pos(all_fragment);
	for(int i=0; i<all_fragment; i++){
		block[i].push_back(i);
		for(int j=0; j<all_fragments[i].Size(); j++){
			block_pos[i].push_back(all_fragments[i].Pos(j));
		}
	}
	
	while(flag2 != 1){
		flag2 = 1;
		
		for(unsigned int i=0; i<block.size(); i++){
			max = -1;
			for(unsigned int j=0; j<block_pos[i].size(); j++){
				if(max < block_pos[i][j]){
					max = block_pos[i][j];
				}
			}
			
			for(unsigned int j=i+1; j<block.size(); j++){
				if(max==-1 || max < block_pos[j][0]){
					continue;
				}
				
				flag = 0;
				for(unsigned int k=0; k<block_pos[j].size(); k++){
					if(find(block_pos[i].begin(), block_pos[i].end(), block_pos[j][k]) != block_pos[i].end()){
						flag = 1;
						flag2 = 0;
						break;
					}
				}
				
				if(flag == 1){
					for(unsigned int k=0; k<block[j].size(); k++){
						block[i].push_back(block[j][k]);
					}
					for(unsigned int k=0; k<block_pos[j].size(); k++){
						if(find(block_pos[i].begin(), block_pos[i].end(), block_pos[j][k]) == block_pos[i].end()){
							block_pos[i].push_back(block_pos[j][k]);
						}
					}
					
					for(unsigned int k=0; k<block_pos[i].size(); k++){
						if(max < block_pos[i][k]){
							max = block_pos[i][k];
						}
					}
					
					block[j].clear();
					block_pos[j].clear();
				}
			}
		}
	}
	
	int block_size = 0;
	vector<int> offset;
	for(unsigned int i=0; i<block.size(); i++){
		if(block[i].size() == 0)	continue;
		sort(block_pos[i].begin(), block_pos[i].end());
		block_size++;
		
		offset.push_back(block_pos[i][0]);
		for(unsigned int j=0; j<block[i].size(); j++){
			all_fragments[block[i][j]].Offset(offset[offset.size()-1]);
		}
	}
	
	hap.resize(block_size);
	for(unsigned int i=0, j=0; i<hap.size(); i++, j++){
		for(; block[j].size()==0; j++)	;
		
		hap[i].Set_fragment(all_fragments, block[j], block_pos[j], offset[i]);
	}
}
