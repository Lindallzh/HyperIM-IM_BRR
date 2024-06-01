//#pragma once
#include <random>
#include <math.h>
//typedef std::mt19937 Myeng; 
//typedef std::binomial_distribution<int, double> Mydist;
//using namespace _gun_cxx;
class HyperGraph
{
private:
    /// _numV: number of nodes in the graph.
    uint32_t _numV;
    /// _numE: number of edges in the graph.
    size_t _numE;
    size_t firstLayerCnt=0;
    size_t totalLayerCnt=0;
    size_t ftRatio=1;
    float ftRatio_1=1;
    /// _numRRsets: number of RR sets.
    size_t _numRRsets = 0;
    std::vector<bool> _vecVisitBool;
    Nodelist _vecVisitNode;
    double pMax=0.0000;
    size_t RR_size = 0;


//typedef std::pair<uint32_t, uint32_t> Layer;

//typedef std::pair<uint32_t, uint32_t> vertexLayer;
/// Edgelist structure from one source/target node
typedef std::vector<uint32_t> layerVertex;
typedef std::vector<uint32_t> ber_position;
/// Graph structure
typedef std::vector<layerVertex> VertexLayer;
typedef std::vector<ber_position> Position_ber;
typedef std::vector<int> neighbor;
typedef std::vector<neighbor> vecNei;
Position_ber _posBer;
VertexLayer _verlay;
vecNei _vecNei;
double pp_possion[8]={0.3679,0.7358,0.9197,0.9810,0.9963,0.9994,0.9999,1.0};
  //  std::vector< std::vector::<std::pair<uint32_t, uint32_t>>> 

    size_t _hit = 0;
    double _numSamplesEval = 0;
    double _hyperedgeAvgEval = 0.0;

    bool _isVanilla = false;

    /// Initialization
    //layer要在这里设置
    void InitHypergraph()
    {
        _numV = (uint32_t)_graph.size();
        //_verlay.resize(_numV);
       // std::cout << " _numV : " << _numV<< std::endl; 
        _verlay.resize(_numV);
        _vecNei.resize(_numV);

      //  for (auto& nbrs : _graph) _numE += nbrs.size();
        
        uint32_t nodeCnt=0;
        
       // std::cout << "InitHypergraph : " << std::endl; 
       // for (auto &nbrs : _graph)
        for (size_t i = 0; i < _graph.size(); i++)

        {
             _numE += _graph[i].size();
             uint32_t ulayerCnt=0;
             uint32_t ulayer=0;
           // std::cout << "vecGRev[i].size() : " <<_graph[i].size()<<std::endl; 

             float preww=0.000000;

            for (size_t j =0; j < _graph[i].size(); j++)
            {
                 
                  _vecNei[i].push_back(j);
                 if(_graph[i][j].second!=preww)
                 {
                    _verlay[i].push_back(ulayerCnt);
                    preww=_graph[i][j].second;
                    if (pMax<_graph[i][j].second&&_graph[i][j].second!=1)
                    {
                        pMax=_graph[i][j].second;
                    }
                
                 }
                 
                 ulayerCnt++; 
                 
               
            }
           // std::random_shuffle(_vecNei[i].begin(),_vecNei[i].end()); 

          // if(ulayerCnt==_graph[i].size()&& _verlay[i].size()==1)
           // {
                //   _verlay[i].push_back(ulayerCnt);
           // }
             
           // std::cout << "nbrs.size() : " <<nbrs.size()<<std::endl; 
        }
       // std::cout << "_numE********************* : " <<_numE<<std::endl;

 for (size_t expand =0; expand < _graph.size(); expand++)
 {
               int totalLayer=_verlay[expand].size();
                int position_1=0;
                int position_2=0;
                 totalLayerCnt= totalLayerCnt+totalLayer;
                // firstLayerCnt=firstLayerCnt+_graph[expand].size();
               if (totalLayer==1)
               {
                   firstLayerCnt=firstLayerCnt+_graph[expand].size();
                  

                   //totalLayerCnt=totalLayerCnt+_graph[expand].size();
               }
               else if (totalLayer>1)
               {
                   firstLayerCnt=firstLayerCnt+_verlay[expand][1];
                  // totalLayerCnt=totalLayerCnt+_graph[expand].size()-_verlay[expand][1];
               }
                
                  for(int ss=1;ss<=totalLayer;ss++)
                  {
                       // std::cout<<"执行到这里了for(int ss=1;ss<=totalLayer;ss++) "<<std::endl;    
                        if (ss==totalLayer)
                          {
                                position_1=_graph[expand].size();
                                position_2=_verlay[expand][ss-1];
                                
                          }
                          else
                          {
                             position_1=_verlay[expand][ss];
                             position_2=_verlay[expand][ss-1];
                          }
          std::random_shuffle(_vecNei[expand].begin()+position_2,_vecNei[expand].begin()+position_1-1);                
                        
      }
      

 }
     //ftRatio=ftRatio*((float)firstLayerCnt/totalLayerCnt);
  ftRatio=(totalLayerCnt/(_graph.size()));//每个顶点的rrset 平均包含几个顶点
 
 ftRatio_1=((float)firstLayerCnt/totalLayerCnt); //每层平均有几个顶点


  // std::cout<<"执行到这里了firstLayerCnt "<<firstLayerCnt<<std::endl;
  // std::cout<<"执行到这里了totalLayerCnt "<<totalLayerCnt<<std::endl;
  //std::cout<<"ftRatio_1 "<<ftRatio_1<<std::endl;

        _FRsets = FRsets(_numV);
        _vecVisitBool = std::vector<bool>(_numV);
        _vecVisitNode = Nodelist(_numV);
        //RefreshHypergraph();
    }

public:
    /// _graph: reverse graph
    const Graph& _graph;
    /// _FRsets: forward cover sets, _FRsets[i] is the node sets that node i can reach
    FRsets _FRsets;
    /// _RRsets: reverse cover sets, _RRsets[i] is the node set that can reach node i
    RRsets _RRsets;

    ProbDist _probDist = WEIGHTS;

    explicit HyperGraph(const Graph& graph) : _graph(graph)
    {
       
        InitHypergraph();
    }

    /// Set cascade model
    void set_prob_dist(const ProbDist dist)
    {
        _probDist = dist;
    }

    void set_vanilla_sample(const bool isVanilla)
    {
        _isVanilla = isVanilla;
    }

    /// Returns the number of nodes in the graph.
    uint32_t get_nodes() const
    {
        return _numV;
    }

    double get_pMax() const
    {
        return pMax;
    }

    /// Returns the number of edges in the graph.
    size_t get_edges() const
    {
        return _numE;
    }

    /// Returns the number of RR sets in the graph.
    size_t get_RR_sets_size() const
    {
        return _numRRsets;
    }
    size_t get_RR_size() const
    {
           return RR_size;
    }
    /// Generate a set of n RR sets
    void BuildRRsetsR2(const size_t numSamples)
    {
        if (numSamples > SIZE_MAX)
        {
            std::cout << "Error:R too large" << std::endl;
            exit(1);
        }

        const auto prevSize = _numRRsets;
        _numRRsets = _numRRsets > numSamples ? _numRRsets : numSamples;
         for (auto i = prevSize; i < numSamples; i++)
            {
                BuildOneRRsetWeightedR2(dsfmt_gv_genrand_uint32_range(_numV), i);
               // BuildOneRRset(dsfmt_gv_genrand_uint32_range(_numV), i);
               //BuildOneRRsetConstant(dsfmt_gv_genrand_uint32_range(_numV), i);
            }
    }
    void BuildRRsets(const size_t numSamples)
    {
        void (*func)(const uint32_t uStart, const size_t hyperIdx);

        if (numSamples > SIZE_MAX)
        {
            std::cout << "Error:R too large" << std::endl;
            exit(1);
        }

        const auto prevSize = _numRRsets;
        _numRRsets = _numRRsets > numSamples ? _numRRsets : numSamples;

        if (_isVanilla)
        {
            // std::cout << "Sample RR set by vanilla method" << std::endl;

            for (auto i = prevSize; i < numSamples; i++)
            {
                BuildOneRRset(dsfmt_gv_genrand_uint32_range(_numV), i);
            }

            return ;
        }

        if (_probDist == WC)
        {
            // std::cout << "Sample RR sets in WC model" << std::endl;

            for (auto i = prevSize; i < numSamples; i++)
            {
                BuildOneRRsetWeighted(dsfmt_gv_genrand_uint32_range(_numV), i);
                 //BuildOneRRset(dsfmt_gv_genrand_uint32_range(_numV), i);
               //BuildOneRRsetConstant(dsfmt_gv_genrand_uint32_range(_numV), i);
            }

            return ;
        }
        else if (_probDist == UNIFORM)
        {
            // std::cout << "Sample RR sets in uniform model" << std::endl;

            for (auto i = prevSize; i < numSamples; i++)
            {
                BuildOneRRsetConstant(dsfmt_gv_genrand_uint32_range(_numV), i);
            }

            return ;
        }
        else if (_probDist == SKEWED || _probDist == WEIGHTS)
        {
            // std::cout << "Sample RR sets in skewed or weights case" << std::endl;

            for (auto i = prevSize; i < numSamples; i++)
            {
                BuildOneRRsetSkewed(dsfmt_gv_genrand_uint32_range(_numV), i);
            }

            return ;
        }
        else
        {
            for (auto i = prevSize; i < numSamples; i++)
            {
                BuildOneRRset(dsfmt_gv_genrand_uint32_range(_numV), i);;
            }
        }

        return ;
    }


    double EvalHyperedgeAvg()
    {
        return _hyperedgeAvgEval;
    }


    void display_hyperedge_stat()
    {
        size_t total_hyperedges = 0;

        for (size_t i = 0; i < _numRRsets; i++)
        {
            total_hyperedges += _RRsets[i].size();
        }

        double ave_hyperedge_size = 1.0 * total_hyperedges / _numRRsets;
        double var = 0.0;

        for (size_t i = 0; i < _numRRsets; i++)
        {
            double diff = _RRsets[i].size() - ave_hyperedge_size;
            var += (diff * diff);
        }

        var = var / _numRRsets;
        std::cout << "final average RR set size: " << ave_hyperedge_size << ", final variance: " << var << std::endl;
        return ;
    }

    double HyperedgeAvg()
    {
        size_t totalHyperedges = 0;
        double avgSize = 0.0;

        for (size_t i = 0; i < _numRRsets; i++)
        {
            totalHyperedges += _RRsets[i].size();
        }

     // return  avgSize = 1.0 * totalHyperedges / _numRRsets;
        return totalHyperedges;
    }

    double HyperedgeMedian()
    {
        std::vector<int> RRsetSize(_numRRsets);
        for (size_t i = 0; i < _numRRsets; i++)
        {
          RRsetSize[i]  = _RRsets[i].size();
          std::cout << RRsetSize[i] << std::endl;
        }

        std::sort(RRsetSize.begin(), RRsetSize.end());
        int index = _numRRsets/2;
        if (_numRRsets%2 == 1)
        {
            return RRsetSize[index];
        }
        else
        {
            return (RRsetSize[index] + RRsetSize[index-1])/2.0;
        }
    } 

    // Generate one RR set
    void BuildOneRRset(const uint32_t uStart, const size_t hyperIdx)
    {
        size_t numVisitNode = 0, currIdx = 0;
        _FRsets[uStart].push_back(hyperIdx);
        _vecVisitNode[numVisitNode++] = uStart;
        _vecVisitBool[uStart] = true;

       // while (currIdx < numVisitNode)
        //{
            const auto expand = _vecVisitNode[currIdx++];

            for (auto& nbr : _graph[expand])
            {   
              
                const auto randDouble = dsfmt_gv_genrand_open_close();
                  const auto nbrId = nbr.first;    
               if (randDouble > nbr.second)
                    continue;
                if (_vecVisitBool[nbrId])
                    continue;               
                     RR_size++;
                     _vecVisitNode[numVisitNode++] = nbrId;
                     _vecVisitBool[nbrId] = true;
                     _FRsets[nbrId].push_back(hyperIdx);
                                
    
            }
       // }
       // std::cout << "RR_size: " <<RR_size << std::endl;

        for (int i = 0; i < numVisitNode; i++) _vecVisitBool[_vecVisitNode[i]] = false;

        _RRsets.push_back(RRset(_vecVisitNode.begin(), _vecVisitNode.begin() + numVisitNode));
    }

    // Independent cascade with weighted probability

   double getPMax(int i)
    {
                int totalLayer=_verlay[i].size();
                int position=_verlay[i][totalLayer-1];
                // The probability int the highest layer
                double p_max= _graph[i][position].second;
                return p_max;
    }
   double getPMin(int i)
    {
        double p_min= _graph[i][0].second;
        return p_min;
    }
  
int factorial(int n) {  // 计算阶乘
    int result = 1;
    for (int i = 1; i <= n; i++) {
        result *= i;
    }
    return result;
}

int combination(int n, int k) {  // 计算组合数
    if (k > n) {
        return 0;
    }
    return factorial(n) / (factorial(k) * factorial(n - k));
}
int position_possion_table(double p_m)
{
 int left = 0, right = 7;
    while (left <= right) {
        int mid = left + (right - left) / 2;
        if (pp_possion[mid] > p_m) {
            right = mid - 1;
        } else {
            left = mid + 1;
        }
    }
    return left;
}
 int creat_possion_table(double p_m, double pos[])
 {
    double p_1,p_2;
    double ff;
    //double pos[];
    double sum;
    p_1=std::pow(1-p_m,4);
    //combination(4,0)
    ff=p_1;
    sum=ff;
    pos[0]=sum;

    p_1=std::pow(1-p_m,3);
    p_2=std::pow(p_m,1); 
    ff=(double)p_1*p_2;
    for(int i=1;i<5;i++)
    {
           sum=sum+ff;
           pos[i]=sum;
    }
     p_1=std::pow(1-p_m,2);
     p_2=std::pow(p_m,2);
     ff=p_1*p_2;

    for(int i=5;i<11;i++)
    {
           sum=sum+ff;
           pos[i+1]=sum;
    }
     p_1=std::pow(1-p_m,1);
     p_2=std::pow(p_m,3);
     ff=p_1*p_2;
    for(int i=11;i<15;i++)
    {
           sum=sum+ff;
           pos[i+1]=sum;
    }
     p_1=std::pow(p_m,4);
     pos[15]=sum+p_1;  
 }

  /*void BuildOneRRsetWeighted(const uint32_t uStart, const size_t hyperIdx)
    {
        size_t numVisitNode = 0, currIdx = 0;
        _FRsets[uStart].push_back(hyperIdx);
        _vecVisitNode[numVisitNode++] = uStart;
        _vecVisitBool[uStart] = true;

        while (currIdx < numVisitNode)
        {
            const auto expand = _vecVisitNode[currIdx++];
            //if (_cascadeModel == IC)
            {
                if (_graph[expand].size() == 0) continue;

                double p =  _graph[expand][0].second;
                double log2Prob = Logarithm(1 - p);

                if (p < 1)
                {
                     int nn=_graph[expand].size();
                     if(nn==1)
                             {
                                const auto nbrId = _graph[expand][0].first; 
                                if (_vecVisitBool[nbrId]) continue; 
                              const auto randDouble = dsfmt_gv_genrand_open_close();
                               if (randDouble>p) continue;
                                                       
                                   _vecVisitNode[numVisitNode++] = nbrId;                        
                                   _vecVisitBool[nbrId] = true;
                                   _FRsets[nbrId].push_back(hyperIdx);
                                

                             }
                            
                           //  if(dsfmt_gv_genrand_open_close()>p_kk)
                          else { //std::cout<<p_m<<std::endl;
                             
                            // if(p_m<1)
                            // {
                             int number=0; 
                            
                            if(nn<=20)
                             {
                               std::default_random_engine generator;
                               std::binomial_distribution<int> distribution(nn,p);
                               number = distribution(generator);
                              
                             }
                             else 
                             {
                               //const auto randDouble = dsfmt_gv_genrand_open_close();
                               //number=position_possion_table(randDouble);
                                 std::default_random_engine generator_1;
                                 std::poisson_distribution<> distribution_possion(1);
                                 number = distribution_possion(generator_1);

                             }                              
                             //随机生成的数有问题，这个要重新生成
                           
                            std::uniform_int_distribution<int> distribution_t(1,nn);
                             
                            for (int ii=0; ii<number; ++ii) {
                                 std::default_random_engine generator_2;
                               int number_1 = distribution_t(generator_2);
                               const auto nbrId = _graph[expand][number_1-1].first; 
                               if (_vecVisitBool[nbrId]) continue;                                     
                                   _vecVisitNode[numVisitNode++] = nbrId;                        
                                   _vecVisitBool[nbrId] = true;
                                   _FRsets[nbrId].push_back(hyperIdx);
                                
                                }
                }
            }
        }
        }

        for (int i = 0; i < numVisitNode; i++) _vecVisitBool[_vecVisitNode[i]] = false;

        _RRsets.push_back(RRset(_vecVisitNode.begin(), _vecVisitNode.begin() + numVisitNode));
    
    }*/

   /*
    void BuildOneRRsetWeighted(const uint32_t uStart, const size_t hyperIdx)
    {
        size_t numVisitNode = 0, currIdx = 0;
        _FRsets[uStart].push_back(hyperIdx);
        _vecVisitNode[numVisitNode++] = uStart;
        _vecVisitBool[uStart] = true;
    
       // while (currIdx < numVisitNode)
        //{
            const auto expand = _vecVisitNode[currIdx++];
            //if (_cascadeModel == IC)
            {
               // if (_graph[expand].size() == 0) continue;

               // double p =  _graph[expand][0].second;
                //

                int totalLayer=_verlay[expand].size();
                int position=_verlay[expand][totalLayer-1];
                // The probability int the highest layer
                double p_m= _graph[expand][position].second;
                //the RR sets must derive from the first layer


               // double log2Prob = Logarithm(1 - p);
               double log2Prob = Logarithm(1 - p_m);

              //  if (p < 1)
               if (p_m < 1)
                {
                    double prob = dsfmt_gv_genrand_open_close();

                    int startPos = Logarithm(prob) / log2Prob;

                    int endPos = _graph[expand].size()-position;

                    while (startPos < endPos)
                    {
                        const auto nbrId = _graph[expand][startPos+position].first;

                        if (_vecVisitBool[nbrId])
                        {
                            int increment = Logarithm(dsfmt_gv_genrand_open_close()) / log2Prob;
                            startPos += (increment + 1);
                            continue;
                        }

                        _vecVisitNode[numVisitNode++] = nbrId;
                        _vecVisitBool[nbrId] = true;
                        _FRsets[nbrId].push_back(hyperIdx);
                        int increment = Logarithm(dsfmt_gv_genrand_open_close()) / log2Prob;
                        startPos += increment + 1;
                        RR_size=RR_size+1;
                    }
                }
                else
                {
                    for (auto& nbr : _graph[expand])
                    {
                        const auto nbrId = nbr.first;

                        if (_vecVisitBool[nbrId])
                            continue;

                        _vecVisitNode[numVisitNode++] = nbrId;
                        _vecVisitBool[nbrId] = true;
                        _FRsets[nbrId].push_back(hyperIdx);
                        RR_size=RR_size+1;
                    }
                }
                for (int i=0;i<_verlay[expand].size()-1;i++)
                {
                      const auto randDouble = dsfmt_gv_genrand_open_close();
                     // int layer=_verlay[expand].size();
                      int position=_verlay[expand][i];
                      int position_re=_verlay[expand][i+1];
                      int ss=position_re-position;
                // The probability int the highest layer
                     double p_k= _graph[expand][position].second;
                     double pp=pow(p_k,ss);
                     if(randDouble>pp)
                     {

                              double log2Prob = Logarithm(1 - p_k);

                        if ( p_k < 1)
                           {
                             double prob = dsfmt_gv_genrand_open_close();

                             int startPos = Logarithm(prob) / log2Prob;
                             
                             int endPos = ss;

                             while (startPos < endPos)
                             {
                                const auto nbrId = _graph[expand][startPos+position].first;
                                if (_vecVisitBool[nbrId])
                                  {
                                     int increment = Logarithm(dsfmt_gv_genrand_open_close()) / log2Prob;
                                     startPos += (increment + 1);
                                     continue;
                                  }

                                _vecVisitNode[numVisitNode++] = nbrId;
                                _vecVisitBool[nbrId] = true;
                                _FRsets[nbrId].push_back(hyperIdx);
                                int increment = Logarithm(dsfmt_gv_genrand_open_close()) / log2Prob;
                                startPos += increment + 1;
                                RR_size=RR_size+1;
                             }
                           }
                           else
                               {
                                   for (auto& nbr : _graph[expand])
                                  {
                                     const auto nbrId = nbr.first;

                                    if (_vecVisitBool[nbrId])
                                      continue;

                                     _vecVisitNode[numVisitNode++] = nbrId;
                                     _vecVisitBool[nbrId] = true;
                                      _FRsets[nbrId].push_back(hyperIdx);
                                      RR_size=RR_size+1;
                                  }
                               }

                            }
                }
            }
       // }

      std::cout << "RR_size: " <<RR_size << std::endl;

        for (int i = 0; i < numVisitNode; i++) _vecVisitBool[_vecVisitNode[i]] = false;

        _RRsets.push_back(RRset(_vecVisitNode.begin(), _vecVisitNode.begin() + numVisitNode));
    }
    */
    
   
    void BuildOneRRsetWeighted(const uint32_t uStart, const size_t hyperIdx)
    {
        size_t numVisitNode = 0, currIdx = 0;
        _FRsets[uStart].push_back(hyperIdx);
        _vecVisitNode[numVisitNode++] = uStart;
        _vecVisitBool[uStart] = true;
            const auto expand = _vecVisitNode[currIdx++];
           int totalLayer=_verlay[expand].size();

             if(totalLayer==1)
                {
                      for (auto& nbr : _graph[expand])
            {
                     const auto nbrId = nbr.first;

                     if (_vecVisitBool[nbrId])
                    continue;

               // const auto randDouble = dsfmt_gv_genrand_open_close();

               // if (randDouble > nbr.second)
                 //   continue;

                _vecVisitNode[numVisitNode++] = nbrId;
                _vecVisitBool[nbrId] = true;
                _FRsets[nbrId].push_back(hyperIdx);
				RR_size++;
            }
                }
        
                if (_graph[expand].size() != 0) 
               { //continue;
                int cntL=(int)( ftRatio_1*(_graph[expand].size()))+1;
                //int totalLayer=_verlay[expand].size();
                int position_1=0;
                int position_2=0;
                  for(int ss=1;ss<=totalLayer;ss++)
                  {
                       // std::cout<<"执行到这里了for(int ss=1;ss<=totalLayer;ss++) "<<std::endl;    
                        if (ss==totalLayer)
                          {
                                position_1=_graph[expand].size();
                                position_2=_verlay[expand][ss-1];
                          }
                          else
                          {
                             position_1=_verlay[expand][ss];
                             position_2=_verlay[expand][ss-1];
                          }
                            // std::cout<<position_1<<" "<<position_2<<std::endl;
                             double p_m= _graph[expand][position_2].second;
                             int nn=position_1-position_2;
                             // double p_kk;
                            // if(nn<=cntL&&ss==1)
                          //  if(nn<=2&&ss<=2)
                            
                             //if(nn<=4&&ss<=2)
                            // if(nn<=4&&ss<=ftRatio)
                             if(nn<=cntL&&ss<=ftRatio)
                             //if(nn==1)
                             {
                              /* const auto nbrId = _graph[expand][position_2].first;
                               if (!_vecVisitBool[nbrId]) 
                               {
                                 _vecVisitNode[numVisitNode++] = nbrId;                        
                                  _vecVisitBool[nbrId] = true;
                                  _FRsets[nbrId].push_back(hyperIdx);

                               }*/
                                for(int kkk=0;kkk<nn;kkk++)
                               {
                                 
                                 
                                 const auto nbrId = _graph[expand][position_2+kkk].first; 
                              if (!_vecVisitBool[nbrId]) 
                                {
                              // const auto randDouble = dsfmt_gv_genrand_open_close();
                             //  if (randDouble>p_m) continue;
                                                       
                                  _vecVisitNode[numVisitNode++] = nbrId;                        
                                  _vecVisitBool[nbrId] = true;
                                  _FRsets[nbrId].push_back(hyperIdx);
                                }
                               } 
                                 //   RR_size=RR_size+1;
                                    //RR_size=RR_size+1;
                              // }
                                

                             }
                            
                           //  if(dsfmt_gv_genrand_open_close()>p_kk)
                          else 
                          { //std::cout<<p_m<<std::endl;
                             
                            // if(p_m<1)
                            // {
                             int number=0; 
                            
                            if(nn<=20)
                             {
                               std::default_random_engine generator;
                               std::binomial_distribution<int> distribution(nn,p_m);
                               number = distribution(generator);
                               //RR_size=RR_size+1;
                               
                             }
                             else 
                             {
                               //const auto randDouble = dsfmt_gv_genrand_open_close();
                               //number=position_possion_table(randDouble);
                                 std::default_random_engine generator_1;
                                 std::poisson_distribution<> distribution_possion(1);
                                 number = distribution_possion(generator_1);
                                 //RR_size=RR_size+1;

                             }                              
                             //随机生成的数有问题，这个要重新生成

                           
                           // std::uniform_int_distribution<int> distribution_t(1,nn);
                           // std::random_shuffle(_vecNei[expand].begin()+position_2,_vecNei[expand].begin()+position_1-1);
                            bool flag=false;
                            for (int ii=0; ii<number; ++ii) {
                              flag=false;
                               //  std::default_random_engine generator_2;
                              // int number_1 = distribution_t(generator_2);
                               //const auto nbrId = _graph[expand][rm[ii]].first; 
                               const auto nbrId = _graph[expand][_vecNei[expand][position_2+ii]].first;
							  // const auto nbrId = _graph[expand][position_1+number_1-1].first; 
                               if (_vecVisitBool[nbrId]) continue;                                     
                                   _vecVisitNode[numVisitNode++] = nbrId;                        
                                   _vecVisitBool[nbrId] = true;
                                   _FRsets[nbrId].push_back(hyperIdx);
                                   flag=true;
                                   //RR_size=RR_size+1;
                                
                                }
                                if (flag)
                                {
                                    RR_size=RR_size+1;
                                }
                  
                  }                    
        }
            
        }
       // }
        // std::cout << "RR_size: " <<RR_size << std::endl;
         // std::cout<<"执行到这里了 "<<std::endl;
        for (int i = 0; i < numVisitNode; i++) _vecVisitBool[_vecVisitNode[i]] = false;
        _RRsets.push_back(RRset(_vecVisitNode.begin(), _vecVisitNode.begin() + numVisitNode));
    }

     void BuildOneRRsetWeightedR2(const uint32_t uStart, const size_t hyperIdx)
    {
        size_t numVisitNode = 0, currIdx = 0;
        _FRsets[uStart].push_back(hyperIdx);
        _vecVisitNode[numVisitNode++] = uStart;
        _vecVisitBool[uStart] = true;
       //buildoneset
       /* while (currIdx<numVisitNode)
        {
           // std::cout<<"currIdx: "<<currIdx<<"numVisitNode: "<<numVisitNode<<" "<<_graph.size()<<std::endl;
            const auto expand = _vecVisitNode[currIdx++];
        
                if (_graph[expand].size() == 0) continue;
                
                int totalLayer=_verlay[expand].size();
               
                if(totalLayer==1)
                {
                      for (auto& nbr : _graph[expand])
            {
                     const auto nbrId = nbr.first;

                     if (_vecVisitBool[nbrId])
                    continue;

                const auto randDouble = dsfmt_gv_genrand_open_close();

                if (randDouble > nbr.second)
                    continue;

                _vecVisitNode[numVisitNode++] = nbrId;
                _vecVisitBool[nbrId] = true;
                _FRsets[nbrId].push_back(hyperIdx);
				RR_size++;
            }
                   
                   
                    continue;

                }

        }*/
                const auto expand = uStart;
                int cntL=(int)(ftRatio*(_graph[expand].size()))+1;
                int totalLayer=_verlay[expand].size();
                int position_1=0;
                int position_2=0;
                for(int ss=1;ss<=totalLayer-1;ss++)
                  {
                       // std::cout<<"执行到这里了for(int ss=1;ss<=totalLayer;ss++) "<<std::endl;    
                        if (ss==totalLayer)
                          {
                                position_1=_graph[expand].size();
                                position_2=_verlay[expand][ss-1];
                          }
                          else
                          {
                             position_1=_verlay[expand][ss];
                             position_2=_verlay[expand][ss-1];
                          }
                            // std::cout<<position_1<<" "<<position_2<<std::endl;
                             double p_m= _graph[expand][position_2].second;
                             int nn=position_1-position_2;
                             // double p_kk;


                             //if(nn<=cntL+2&&ss<=3)
                            // if(ss<=2)
                             if(nn>=cntL&&ss<=ftRatio)
                            //if(ss>=2&&nn==2)
                            // if(nn<=cntL+2)
                            // if(nn==1)
                             {
                              /* const auto nbrId = _graph[expand][position_2].first;
                               if (!_vecVisitBool[nbrId]) 
                               {
                                 _vecVisitNode[numVisitNode++] = nbrId;                        
                                  _vecVisitBool[nbrId] = true;
                                  _FRsets[nbrId].push_back(hyperIdx);

                               }*/
                                for(int kkk=0;kkk<nn;kkk++)
                               {
                                 
                                 
                                 const auto nbrId = _graph[expand][position_2+kkk].first;
                               //  const auto randDouble = dsfmt_gv_genrand_open_close();
                             // if (randDouble>p_m) continue; 
                              if (!_vecVisitBool[nbrId]) 
                                {
                               
                                                       
                                  _vecVisitNode[numVisitNode++] = nbrId;                        
                                  _vecVisitBool[nbrId] = true;
                                  _FRsets[nbrId].push_back(hyperIdx);
                                }
                               } 
                                 //   RR_size=RR_size+1;
                                    //RR_size=RR_size+1;
                              // }
                                

                             }
                            


/*
                             if(nn==1)
                             {
                                const auto nbrId = _graph[expand][position_2].first; 
                                if (_vecVisitBool[nbrId]) continue; 
                              //  const auto randDouble = dsfmt_gv_genrand_open_close();
                             //  if (randDouble>p_m) continue;
                                                       
                                   _vecVisitNode[numVisitNode++] = nbrId;                        
                                   _vecVisitBool[nbrId] = true;
                                   _FRsets[nbrId].push_back(hyperIdx);
                                

                             }*/
                            
                           //  if(dsfmt_gv_genrand_open_close()>p_kk)
                          else { //std::cout<<p_m<<std::endl;
                             
                            // if(p_m<1)
                            // {
                             int number=0; 
                            
                            if(nn<=20)
                             {
                               std::default_random_engine generator;
                               std::binomial_distribution<int> distribution(nn,p_m);
                               number = distribution(generator);
							   RR_size=RR_size+1;
                              
                             }
                             else 
                             {
                               //const auto randDouble = dsfmt_gv_genrand_open_close();
                               //number=position_possion_table(randDouble);
                                 std::default_random_engine generator_1;
                                 std::poisson_distribution<> distribution_possion(1);
                                 number = distribution_possion(generator_1);
								 RR_size=RR_size+1;

                             }                              
                             //随机生成的数有问题，这个要重新生成

                           
                           // std::uniform_int_distribution<int> distribution_t(1,nn);
                           // std::random_shuffle(_vecNei[expand].begin()+position_2,_vecNei[expand].begin()+position_1-1);
                            for (int ii=0; ii<number; ++ii) {
                               //  std::default_random_engine generator_2;
                              // int number_1 = distribution_t(generator_2);
                               //const auto nbrId = _graph[expand][rm[ii]].first; 
                               const auto nbrId = _graph[expand][_vecNei[expand][position_2+ii]].first;
							  // const auto nbrId = _graph[expand][position_1+number_1-1].first; 
                               if (_vecVisitBool[nbrId]) continue;                                     
                                   _vecVisitNode[numVisitNode++] = nbrId;                        
                                   _vecVisitBool[nbrId] = true;
                                   _FRsets[nbrId].push_back(hyperIdx);
                                
                                }
                  
                  }                    
        
    }
        //}
         // std::cout<<"执行到这里了 "<<std::endl;
        for (int i = 0; i < numVisitNode; i++) _vecVisitBool[_vecVisitNode[i]] = false;
        _RRsets.push_back(RRset(_vecVisitNode.begin(), _vecVisitNode.begin() + numVisitNode));
    }

    /* independent cascade with constant probability */
    void BuildOneRRsetConstant(const uint32_t uStart, const size_t hyperIdx)
    {
        size_t numVisitNode = 0, currIdx = 0;
        _FRsets[uStart].push_back(hyperIdx);
        _vecVisitNode[numVisitNode++] = uStart;
        _vecVisitBool[uStart] = true;
        const double p =  _graph[0][0].second;
        const double const_prob = Logarithm(1 - p);

        if (p == 1)
        {
            while (currIdx < numVisitNode)
            {
                const auto expand = _vecVisitNode[currIdx++];

                //std::cout<<_graph[expand].size()<<std::endl;
                if (_graph[expand].size() == 0) continue;

                for (auto& nbr : _graph[expand])
                {
                    const auto nbrId = nbr.first;

                    //std::cout<<nbr.first<<" "<<nbr.second<<std::endl;
                    if (_vecVisitBool[nbrId])
                        continue;

                    _vecVisitNode[numVisitNode++] = nbrId;
                    _vecVisitBool[nbrId] = true;
                    _FRsets[nbrId].push_back(hyperIdx);
                }
            }
        }
        else
        {
            while (currIdx < numVisitNode)
            {
                const auto expand = _vecVisitNode[currIdx++];

                if (0 == _graph[expand].size())
                {
                    continue;
                }

                int startPos = Logarithm(dsfmt_gv_genrand_open_close()) / const_prob;
                int endPos = _graph[expand].size();

                while (startPos < endPos)
                {
                    //std::cout<<"enter loop"<<std::endl;
                    const auto nbrId = _graph[expand][startPos].first;

                    if (!_vecVisitBool[nbrId])
                    {
                        _vecVisitNode[numVisitNode++] = nbrId;
                        _vecVisitBool[nbrId] = true;
                        _FRsets[nbrId].push_back(hyperIdx);
                    }

                    int increment = Logarithm(dsfmt_gv_genrand_open_close()) / const_prob;
                    startPos += increment + 1;
                }
            }
        }

        for (int i = 0; i < numVisitNode; i++) _vecVisitBool[_vecVisitNode[i]] = false;

        _RRsets.push_back(RRset(_vecVisitNode.begin(), _vecVisitNode.begin() + numVisitNode));
    }

    // independent cascade with skewed distribution
    // independent cascade with skewed distribution 
    void BuildOneRRsetSkewed(const uint32_t uStart, const size_t hyperIdx)
    {
        size_t numVisitNode = 0, currIdx = 0;
        _FRsets[uStart].push_back(hyperIdx);
        _vecVisitNode[numVisitNode++] = uStart;
        _vecVisitBool[uStart] = true;
        double p_threshold = 0.1;

        while (currIdx < numVisitNode)
        {
            const auto expand = _vecVisitNode[currIdx++];
            size_t out_degree = _graph[expand].size();

            if (out_degree > 0)
            {
                size_t startMin = 0;
                size_t endMax = out_degree;

                while (startMin < endMax)
                {
                    const auto &currentedge = _graph[expand][startMin];
                    const auto node_prob = currentedge.second;
                    const auto nbrId = currentedge.first;

                    if (node_prob < p_threshold)
                    {
                        break;
                    }

                    const auto randDouble = dsfmt_gv_genrand_open_close();
                    startMin++;

                    if (randDouble > node_prob) continue;

                    if (_vecVisitBool[nbrId]) continue;

                    _vecVisitNode[numVisitNode++] =  nbrId;
                    _vecVisitBool[nbrId]  = true;
                    _FRsets[nbrId].push_back(hyperIdx);
                }

                while (startMin < endMax)
                {
                    double bucket_probability = _graph[expand][startMin].second;
                    const double log_prob = Logarithm(1 - bucket_probability);
                    double prob = dsfmt_gv_genrand_open_close();
                    startMin += floor(Logarithm(prob) / log_prob);

                    if (startMin >= endMax)
                    {
                        break;
                    }

                    const auto &currentedge = _graph[expand][startMin];
                    const auto nbrId = currentedge.first;
                    const auto accept_probability = currentedge.second;
                    double randDouble = dsfmt_gv_genrand_open_close();
                    startMin++;

                    if (randDouble > accept_probability / bucket_probability || _vecVisitBool[nbrId])
                    {
                        continue;
                    }

                    _vecVisitNode[numVisitNode++] = nbrId;
                    _vecVisitBool[nbrId] = true;
                    _FRsets[nbrId].push_back(hyperIdx);
                }
            }
        }

        for (int i = 0; i < numVisitNode; i++) _vecVisitBool[_vecVisitNode[i]] = false;

        _RRsets.push_back(RRset(_vecVisitNode.begin(), _vecVisitNode.begin() + numVisitNode));
    }

    // Evaluate the influence spread of a seed set on current generated RR sets
    double CalculateInf(const Nodelist& vecSeed)
    {
        std::vector<bool> vecBoolVst = std::vector<bool>(_numRRsets);
       // std::vector<bool> vecBoolSeed(_numV);

       // for (auto seed : vecSeed) vecBoolSeed[seed] = true;

        for (auto seed : vecSeed)
        {
            for (auto node : _FRsets[seed])
            {
                vecBoolVst[node] = true;
            }
        }

        size_t count = std::count(vecBoolVst.begin(), vecBoolVst.end(), true);
        return 1.0 * count * _numV / _numRRsets;
    }

    // Efficiently estimate the influence spread with sampling error epsilon within probability 1-delta
    double EfficInfVldtAlg(const Nodelist& vecSeed, const double delta = 1e-3, const double eps = 0.01)
    {
        const double c = 2.0 * (exp(1.0) - 2.0);
        const double LambdaL = 1.0 + 2.0 * c * (1.0 + eps) * log(2.0 / delta) / (eps * eps);
        size_t numHyperEdge = 0;
        size_t numCoverd = 0;
        std::vector<bool> vecBoolSeed(_numV);

        for (auto seed : vecSeed) vecBoolSeed[seed] = true;

        while (numCoverd < LambdaL)
        {
            numHyperEdge++;
            size_t numVisitNode = 0, currIdx = 0;
            const auto uStart = dsfmt_gv_genrand_uint32_range(_numV);

            if (vecBoolSeed[uStart])
            {
                // Stop, this sample is covered
                numCoverd++;
                continue;
            }

            _vecVisitNode[numVisitNode++] = uStart;
            _vecVisitBool[uStart] = true;

            while (currIdx < numVisitNode)
            {
                const auto expand = _vecVisitNode[currIdx++];

                for (auto& nbr : _graph[expand])
                {
                    const auto nbrId = nbr.first;

                    if (_vecVisitBool[nbrId])
                        continue;

                    const auto randDouble = dsfmt_gv_genrand_open_close();

                    if (randDouble > nbr.second)
                        continue;

                    if (vecBoolSeed[nbrId])
                    {
                        // Stop, this sample is covered
                        numCoverd++;
                        goto postProcess;
                    }

                    _vecVisitNode[numVisitNode++] = nbrId;
                    _vecVisitBool[nbrId] = true;
                }
            }

postProcess:

            for (auto i = 0; i < numVisitNode; i++)
                _vecVisitBool[_vecVisitNode[i]] = false;
        }

        return 1.0 * numCoverd * _numV / numHyperEdge;
    }

    // Refresh the hypergraph
    void RefreshHypergraph()
    {
        if (_RRsets.size() != 0)
        {
            for (auto i = _numRRsets; i--;)
            {
                RRset().swap(_RRsets[i]);
            }

            RRsets().swap(_RRsets);

            for (auto i = _numV; i--;)
            {
                FRset().swap(_FRsets[i]);
            }
        }

        _numRRsets = 0;
        _hit = 0;
    }

    // Release memory
    void ReleaseMemory()
    {
        RefreshHypergraph();
        std::vector<bool>().swap(_vecVisitBool);
        Nodelist().swap(_vecVisitNode);
        FRsets().swap(_FRsets);
    }

    void BuildOneRRsetEarlyStopByVanilla(std::unordered_set<uint32_t> &connSet, const uint32_t uStart, const size_t hyperIdx)
    {
        size_t numVisitNode = 0, currIdx = 0;
        _FRsets[uStart].push_back(hyperIdx);
        _vecVisitNode[numVisitNode++] = uStart;
        _vecVisitBool[uStart] = true;

        if (connSet.find(uStart) != connSet.end())
        {
            _hit++;
            goto finished;
        }

        while (currIdx < numVisitNode)
        {
            const auto expand = _vecVisitNode[currIdx++];

            if (0 == _graph[expand].size())
            {
                continue;
            }

            for (auto& nbr : _graph[expand])
            {
                const auto nbrId = nbr.first;

                if (_vecVisitBool[nbrId])
                    continue;

                const auto randDouble = dsfmt_gv_genrand_open_close();

                if (randDouble > nbr.second)
                    continue;

                _vecVisitNode[numVisitNode++] = nbrId;
                _vecVisitBool[nbrId] = true;
                _FRsets[nbrId].push_back(hyperIdx);
                RR_size++;

                if (connSet.find(nbrId) != connSet.end())
                {
                    _hit++;
                    goto finished;
                }
            }
        }



finished:

        for (int i = 0; i < numVisitNode; i++) _vecVisitBool[_vecVisitNode[i]] = false;

        _RRsets.push_back(RRset(_vecVisitNode.begin(), _vecVisitNode.begin() + numVisitNode));
    }


    void BuildOneRRsetEarlyStopBySubsim(std::unordered_set<uint32_t> &connSet, const uint32_t uStart, const size_t hyperIdx)
    {
        size_t numVisitNode = 0, currIdx = 0;
        _FRsets[uStart].push_back(hyperIdx);
        _vecVisitNode[numVisitNode++] = uStart;
        _vecVisitBool[uStart] = true;

        if (connSet.find(uStart) != connSet.end())
        {
            _hit++;
            goto finished;
        }

        while (currIdx < numVisitNode)
        {
            const auto expand = _vecVisitNode[currIdx++];
            const double p =  _graph[expand][0].second;

            if (0 == _graph[expand].size())
            {
                continue;
            }

            if (p >= 1.0)
            {
                for (auto& nbr : _graph[expand])
                {
                    const auto nbrId = nbr.first;

                    if (_vecVisitBool[nbrId])
                        continue;

                    _vecVisitNode[numVisitNode++] = nbrId;
                    _vecVisitBool[nbrId] = true;
                    _FRsets[nbrId].push_back(hyperIdx);
                    RR_size++;
                }

                continue;
            }

            const double const_prob = Logarithm(1 - p);
            int startPos = Logarithm(dsfmt_gv_genrand_open_close()) / const_prob;
            int endPos = _graph[expand].size();

            while (startPos < endPos)
            {
                const auto nbrId = _graph[expand][startPos].first;

                if (!_vecVisitBool[nbrId])
                {
                    _vecVisitNode[numVisitNode++] = nbrId;
                    _vecVisitBool[nbrId] = true;
                    _FRsets[nbrId].push_back(hyperIdx);
                    RR_size++;
                }

                if (connSet.find(nbrId) != connSet.end())
                {
                    _hit++;
                    goto finished;
                }

                int increment = Logarithm(dsfmt_gv_genrand_open_close()) / const_prob;
                startPos += increment + 1;
            }
        }

finished:

        for (int i = 0; i < numVisitNode; i++) _vecVisitBool[_vecVisitNode[i]] = false;

        _RRsets.push_back(RRset(_vecVisitNode.begin(), _vecVisitNode.begin() + numVisitNode));
    }

    void BuildRRsetsEarlyStop(std::unordered_set<uint32_t> &connSet, const int numSamples)
    {
        const auto prevSize = _numRRsets;
        _numRRsets = _numRRsets > numSamples ? _numRRsets : numSamples;

        if (_isVanilla)
        {
            // std::cout << "Sample RR sets with early stop by Vanilla method" << std::endl;
            for (auto i = prevSize; i < numSamples; i++)
            {
                BuildOneRRsetEarlyStopByVanilla(connSet, dsfmt_gv_genrand_uint32_range(_numV), i);
            }
        }
        else
        {
            // std::cout << "Sample RR sets with early stop By SUBSIM" << std::endl;
            for (auto i = prevSize; i < numSamples; i++)
            {
                BuildOneRRsetEarlyStopBySubsim(connSet, dsfmt_gv_genrand_uint32_range(_numV), i);
            }
        }
    }

    double EvalSeedSetInfByVanilla(std::unordered_set<uint32_t> &connSet, const int numSamples)
    {
        uint32_t numCovered = 0;
        int64_t totalHyperedgeSize = 0;
        _numSamplesEval = numSamples;

        for (int i = 1; i < numSamples; i++)
        {
            uint32_t uStart = dsfmt_gv_genrand_uint32_range(_numV);
            size_t numVisitNode = 0, currIdx = 0;
            _vecVisitNode[numVisitNode++] = uStart;
            _vecVisitBool[uStart] = true;

            if (connSet.find(uStart) != connSet.end())
            {
                numCovered++;
                goto finished;
            }

            while (currIdx < numVisitNode)
            {
                const auto expand = _vecVisitNode[currIdx++];

                if (0 == _graph[expand].size())
                {
                    continue;
                }

                const double p =  _graph[expand][0].second;

                for (auto& nbr : _graph[expand])
                {
                    const auto nbrId = nbr.first;

                    if (_vecVisitBool[nbrId])
                        continue;

                    const auto randDouble = dsfmt_gv_genrand_open_close();

                    if (randDouble > nbr.second)
                        continue;

                    _vecVisitNode[numVisitNode++] = nbrId;
                    _vecVisitBool[nbrId] = true;
                    RR_size++;

                    if (connSet.find(nbrId) != connSet.end())
                    {
                        numCovered++;
                        goto finished;
                    }
                }
            }

finished:
            totalHyperedgeSize += numVisitNode;

            for (int i = 0; i < numVisitNode; i++) _vecVisitBool[_vecVisitNode[i]] = false;
        }

        _hyperedgeAvgEval = 1.0 * totalHyperedgeSize / numSamples;
        return 1.0 * numCovered * _numV / numSamples;
    }

    double EvalSeedSetInfBySubsim(std::unordered_set<uint32_t> &connSet, const int numSamples)
    {
        uint32_t numCovered = 0;
        int64_t totalHyperedgeSize = 0;
        _numSamplesEval = numSamples;

        for (int i = 1; i < numSamples; i++)
        {
            uint32_t uStart = dsfmt_gv_genrand_uint32_range(_numV);
            size_t numVisitNode = 0, currIdx = 0;
            _vecVisitNode[numVisitNode++] = uStart;
            _vecVisitBool[uStart] = true;

            if (connSet.find(uStart) != connSet.end())
            {
                numCovered++;
                goto finished;
            }

            while (currIdx < numVisitNode)
            {
                const auto expand = _vecVisitNode[currIdx++];

                if (0 == _graph[expand].size())
                {
                    continue;
                }

                const double p =  _graph[expand][0].second;

                if (p >= 1.0)
                {
                    for (auto& nbr : _graph[expand])
                    {
                        const auto nbrId = nbr.first;

                        if (_vecVisitBool[nbrId])
                            continue;

                        _vecVisitNode[numVisitNode++] = nbrId;
                        _vecVisitBool[nbrId] = true;

                        if (connSet.find(nbrId) != connSet.end())
                        {
                            numCovered++;
                            goto finished;
                        }
                    }

                    continue;
                }

                const double const_prob = Logarithm(1 - p);
                int startPos = Logarithm(dsfmt_gv_genrand_open_close()) / const_prob;
                int endPos = _graph[expand].size();

                while (startPos < endPos)
                {
                    const auto nbrId = _graph[expand][startPos].first;

                    if (!_vecVisitBool[nbrId])
                    {
                        _vecVisitNode[numVisitNode++] = nbrId;
                        _vecVisitBool[nbrId] = true;
                        RR_size++;
                    }

                    if (connSet.find(nbrId) != connSet.end())
                    {
                        numCovered++;
                        goto finished;
                    }

                    int increment = Logarithm(dsfmt_gv_genrand_open_close()) / const_prob;
                    startPos += increment + 1;
                }
            }

finished:
            totalHyperedgeSize += numVisitNode;

            for (int i = 0; i < numVisitNode; i++) _vecVisitBool[_vecVisitNode[i]] = false;
        }

        _hyperedgeAvgEval = 1.0 * totalHyperedgeSize / numSamples;
        return 1.0 * numCovered * _numV / numSamples;
    }

    double EvalSeedSetInf(std::unordered_set<uint32_t> &connSet, const int numSamples)
    {
        if (_isVanilla)
        {
            return EvalSeedSetInfByVanilla(connSet, numSamples);
        }
        else
        {
            return EvalSeedSetInfBySubsim(connSet, numSamples);
        }
    }

    double CalculateInfEarlyStop()
    {
        return 1.0 * _hit * _numV / _numRRsets;
    }
};

using THyperGraph = HyperGraph;
using PHyperGraph = std::shared_ptr<THyperGraph>;
