#pragma once
#include <map>
#include <vector>
class GraphBase
{
public:
    /// Format the input for future computing, which is much faster for loading. Vector serialization is used.
   static bool cmp2(Edge a, Edge b)
{
     return a.second<b.second;
}
   
static bool greater_first(const std::pair<uint32_t, float> x, const std::pair<uint32_t, float> y)
{
    if (x.second > y.second) return true;
    else if (x.second < y.second) return false;
    else if (x.first < y.first) return true;
    else return false;
}
    static void FormatGraph(const std::string filename, ProbDist probDist, const float sum, const float prob, const std::string skewType)
    {
        size_t numV, numE;
        uint32_t srcId, dstId;
        float weight = 0;
        size_t firstlayer=0;
       // std::cout << "The file \"" + filename + "\" is opened\n";
        std::ifstream infile(filename);

        if (!infile.is_open())
        {
            std::cout << "The file \"" + filename + "\" can NOT be opened\n";
            return;
        }
      // std::cout << "##########Formate the graph**********" << std::endl;
        infile >> numV >> numE;
        Graph vecGRev(numV);

        //Graph_Nei vecNei(numV);
       // vertexLayer verLayer(numV);
        std::vector<size_t> vecInDeg(numV);

        for (auto i = numE; i--;)
        {
            //if (probDist == WEIGHTS)
            if (probDist == WC)
            {
                infile >> srcId >> dstId >> weight;
            }
            else
            {
                infile >> srcId >> dstId;
            }

            vecGRev[dstId].push_back(Edge(srcId, weight));
           
        }

        infile.close();

        for (auto idx = 0; idx < numV; idx++)
        {
            vecInDeg[idx] = vecGRev[idx].size();
        }

       /* if (probDist == WC)
        {
            for (size_t i = 0; i < vecGRev.size(); i++)
            {
                if (vecGRev[i].size() == 0) continue;

                weight = sum / vecInDeg[i];

                for (size_t j = 0; j < vecGRev[i].size(); j++)
                {
                    vecGRev[i][j].second = weight;
                }
            }
        }*/
        // sort(vecGRev[i][j].begin(), vecGRev[i][j].end());
        //重新设计分层算法
      
      //std::cout << "The file \"" + filename + "\" can NOT be opened\n";
      //std::cout << "Formate the graph**********" << std::endl;

    
      // float prop=0.00000;
      // std::poisson_distribution<> d(1);
      // 分组
       if (probDist == WC)
        {
             // std::cout << "Formate the graph**********" << std::endl;
             // nodeGroup->clear();
             // nodeCnt.clear();
            std::map<float, int> *nodeGroup=new std::map<float,int>;
            std::vector<float> nodeCnt;
             for (size_t i = 0; i < vecGRev.size(); i++)
            {
                if (vecGRev[i].size() == 0) continue;
                
               // std::cout << "vecGRev[i].size() " <<vecGRev[i].size()<<std::endl;

                sort(vecGRev[i].begin(),vecGRev[i].end(),cmp2);
               // std::cout << "vecGRev[i].size()++++ " <<vecGRev[i].size()<<std::endl;
                 for (size_t j = 0; j < vecGRev[i].size(); j++)
                 {
                    (*nodeGroup)[vecGRev[i][j].second]++;
                 }
                // 
                float ii = 0.0;
                size_t j=1;
                size_t totalCnt=0;
                size_t front=0;
                
               // verLayer[i]=(*nodeGroup).size();
               //这是分组情况
                for(std::map<float,int>::iterator iter=(*nodeGroup).begin(); iter != (*nodeGroup).end() ; ++iter)
                {
                        nodeCnt.push_back(ii+(float)(iter->second));
                       // [i][j]
                       // verLayer[i].push_back(std::make_pair(j, front));
                        //ii=ii+1.0; 
                        totalCnt=totalCnt+iter->second;
                        //j++;
                        
                }
              //  std::cout << "vecGRev[i].size()++++333333 " <<vecGRev[i].size()<<std::endl;
                 std::sort(nodeCnt.begin(),nodeCnt.end());
       //  std::cout << "vecGRev[i].size()++++ 2" <<vecGRev[i].size()<<std::endl;
                 size_t t = 0;
                 size_t tt=0;

               
                 for(std::map<float,int>::iterator iter=(*nodeGroup).begin(); iter != (*nodeGroup).end() ; ++iter)
                 {
                      //= nodeCnt[t] / totalCnt;
                    for(size_t k = 0; k < iter->second; k++)
                    {
                           vecGRev[i][tt].second = nodeCnt[t] / (totalCnt*iter->second); 
                           tt++; 
                          // k++;
                    }
                    t++;
                    
                 }
                  sort(vecGRev[i].begin(), vecGRev[i].end(), greater_first);
                 (*nodeGroup).clear();
                 nodeCnt.clear();
            }

        }
        /*if (probDist == WC)
        {
             // std::cout << "Formate the graph**********" << std::endl;
             // nodeGroup->clear();
             // nodeCnt.clear();
            std::map<float, int> *nodeGroup=new std::map<float,int>;
            std::vector<float> nodeCnt;
             for (size_t i = 0; i < vecGRev.size(); i++)
            {
                if (vecGRev[i].size() == 0) continue;
               
               // std::cout << "vecGRev[i].size() " <<vecGRev[i].size()<<std::endl;

                sort(vecGRev[i].begin(),vecGRev[i].end(),cmp2);
                

               // std::cout << "vecGRev[i].size()++++ " <<vecGRev[i].size()<<std::endl;
                 for (size_t j = 0; j < vecGRev[i].size(); j++)
                 {
                    (*nodeGroup)[vecGRev[i][j].second]++;
                 }
                // 
                float ii = 1.0;
                size_t j=1;
                size_t totalCnt=0;
                size_t front=0;
                
               // verLayer[i]=(*nodeGroup).size();
               //这是分组情况
                for(std::map<float,int>::iterator iter=(*nodeGroup).begin(); iter != (*nodeGroup).end() ; ++iter)
                {
                        nodeCnt.push_back(ii+(float)(iter->second));
                       // [i][j]
                       // verLayer[i].push_back(std::make_pair(j, front));
                        ii=ii+1.0; 
                        totalCnt=totalCnt+j+iter->second;
                        j++;
                        
                }
              //  std::cout << "vecGRev[i].size()++++333333 " <<vecGRev[i].size()<<std::endl;
                 std::sort(nodeCnt.begin(),nodeCnt.end());
       //  std::cout << "vecGRev[i].size()++++ 2" <<vecGRev[i].size()<<std::endl;
                 size_t t = 0;
                 size_t tt=0;
				 size_t ll=1;

               float frac_1=(float) (1.000)/ (log2((*nodeGroup).size())+1);
                 for(std::map<float,int>::iterator iter=(*nodeGroup).begin(); iter != (*nodeGroup).end() ; ++iter)
                 {
                      //= nodeCnt[t] / totalCnt;
                      float frac_2=(float) 1.000/iter->second;
                    for(size_t k = 0; k < iter->second; k++)
                    {
                          // vecGRev[i][tt].second = nodeCnt[t] / (totalCnt*iter->second);
						   float frac_3=(float) 1.000/ll;
                           // vecGRev[i][k].second = frac_1*frac_2*frac_3; 
                            vecGRev[i][k].second = frac_1*frac_2; 
                           tt++; 
                          // k++;
                    }
					ll++;
                   // t++;
                    
                 }
                  sort(vecGRev[i].begin(), vecGRev[i].end(), greater_first);
                 (*nodeGroup).clear();
                 nodeCnt.clear();
            }

        }*/
        else if (probDist == UNIFORM)
        {
            // Uniform probability
            for (auto& nbrs : vecGRev)
            {
                for (auto& nbr : nbrs)
                {
                    nbr.second = prob;
                }
            }
        }
        // exponential distribution with lambada equal to its in-degree
        else if (probDist == SKEWED)
        {

            if (skewType == "weibull")
            {
                std::default_random_engine generator(time(NULL));
                double min_value = (1e-8 < 1.0/numV)? 1e-8: 1.0/numV;

                for (size_t i = 0; i< vecGRev.size(); i++)
                {

                    if (vecGRev[i].size() == 0) continue;
                    double sum = 0.0;
                    for (size_t j = 0; j < vecGRev[i].size(); j++)
                    {
                        // random number from (0, 10)
                        double a = dsfmt_gv_genrand_open_open() * 10;
                        double b = dsfmt_gv_genrand_open_open() * 10;
                        std::weibull_distribution<double> distribution(a,b);
                        auto weight =  distribution(generator);
                        vecGRev[i][j].second = weight;

                        sum += weight;
                    }

                    for (size_t j = 0; j < vecGRev[i].size(); j++)
                    {
                        auto weight = vecGRev[i][j].second/sum;
                        vecGRev[i][j].second = (weight > min_value)?weight:min_value; 
                    }
                    sort(vecGRev[i].begin(), vecGRev[i].end(), greater_first);
                }
            }
            else
            {
                double min_value = (1e-8 < 1.0/numV)? 1e-8: 1.0/numV;
                for (size_t i = 0; i<vecGRev.size(); i++)
                {
                    if (vecGRev[i].size() == 0) continue;
                    double sum = 0.0;
                    for (size_t j = 0; j < vecGRev[i].size(); j++)
                    {
                        // lambda = 1
                        auto weight = -log ( 1.0 - dsfmt_gv_genrand_open_open());
                        vecGRev[i][j].second = weight;
                        sum += weight;
                    }

                    for (size_t j = 0; j < vecGRev[i].size(); j++)
                    {
                        double weight = vecGRev[i][j].second/sum;
                        vecGRev[i][j].second = (weight > min_value)?weight:min_value;
                    }
                    sort(vecGRev[i].begin(), vecGRev[i].end(), greater_first);
                }
            }
        }
       /* else if (probDist == WEIGHTS)
        {
            for (size_t i = 0; i < vecGRev.size(); i++)
            {
                sort(vecGRev[i].begin(), vecGRev[i].end(), greater_first);
            }
        }*/

       /* for (size_t i = 0; i < vecGRev.size(); i++)
        {
             for (size_t j = 0; j < vecGRev[i].size(); j++)
             {
                 std::cout << "i " << i<<" "<<vecGRev[i][j].first<<" "<<vecGRev[i][j].second<<std::endl;
             }
        }*/

        std::cout << "probability distribution: " << probDist << std::endl;
        TIO::SaveGraphStruct(filename, vecGRev, true);
        TIO::SaveGraphProbDist(filename, (int)probDist);
        std::cout << "The graph is formatted!" << std::endl;
    }

    /// Load graph via vector deserialization.
    static void LoadGraph(Graph &graph, const std::string graphName)
    {
        TIO::LoadGraphStruct(graphName, graph, true);
        return ;
    }

    static int LoadGraphProbDist(const std::string graphName)
    {
        return TIO::LoadGraphProbDist(graphName);
    }
};
