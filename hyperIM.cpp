#include "stdafx.h"
#include "SFMT/dSFMT/dSFMT.c"
#include "alg.cpp"

void init_random_seed()
{
    // Randomize the seed for generating random numbers
    dsfmt_gv_init_gen_rand(static_cast<uint32_t>(time(nullptr)));
}

int main(int argc, char* argv[])
{

    TArgument Arg(argc, argv);

    if (Arg._probDist == PROB_DIST_ERROR)
    {
        LogInfo("The input probability distribution is not supported:", Arg._probDistStr);
        LogInfo("The supported probability distribution: weights, wc, uniform, skewed");
        return 0;
    }

    if (Arg._func == FUNC_ERROR)
    {
        LogInfo("The input func is not supported: ", Arg._funcStr);
        LogInfo("The supported func: format, im");
    }

    init_random_seed();

    std:: string dataname;
    std:: string inputpath;

     for(int i=1; i<argc ; i++){
        string input = argv[i];
        if(input.compare("--inputpath") == 0) inputpath = argv[++i];
        else if(input.compare("--dataname") == 0) dataname = argv[++i];
     }
    // the input of the hypergraph 
    HyperGraphLoad *graphLoad = new HyperGraphLoad(inputpath, dataname);
   // format the graph: setting the weights
    //const std::string infilename = Arg._dir + "/" + Arg._graphname;
    const std::string infilename = Arg._dir + "/" + dataname;
    // path = inputpath + dataname + ".txt";
    if (Arg._func == FORMAT)
    {
        // Format the graph
       //  std::cout << "Formate the graph" << std::endl;
        GraphBase::FormatGraph(infilename, Arg._probDist, Arg._wcVar, Arg._probEdge, Arg._skewType);
       //  GraphBase::FormatGraph(dataname, Arg._probDist, Arg._wcVar, Arg._probEdge, Arg._skewType);
        return 0;
    }

    Timer mainTimer("main");
    // Load the reverse graph
    Graph graph;
  //  std::cout << "Load the Graph" << std::endl;
  // 在load里花费时间比较长
    GraphBase::LoadGraph(graph, infilename);
   // GraphBase::LoadGraph(graph, dataname);
    int probDist = GraphBase::LoadGraphProbDist(infilename);
   //  std::cout << "finish Load the Graph++++" << std::endl;
   //int probDist = GraphBase::LoadGraphProbDist(dataname);

    // Initialize a result object to record the results
    TResult tRes;
   // std::cout << "finish Load the Graph6666666666" << std::endl;
    TAlg tAlg(graph, tRes);
  //  std::cout << "finish Load the Graph++++&&&&&&&&&&&" << std::endl;
    tAlg.set_vanilla_sample(Arg._vanilla);
    tAlg.set_prob_dist((ProbDist)probDist); // Set propagation model
    auto delta = Arg._delta;

    if (delta < 0) delta = 1.0 / graph.size();
  
    int seedSize = Arg._seedsize;
    std::cout << "seedSize k=" << seedSize << std::endl;
    Arg.build_outfilename(seedSize, (ProbDist)probDist, graph);
    std::cout << "---The Begin of " << Arg._outFileName << "---\n";
    time_t start = 0,end = 0;  
    time(&start);
    if (!Arg._hist)
    {
        tAlg.subsimOnly(seedSize, Arg._eps, delta);
    }
    else
    {
        std::cout <<"HIST is invoked." <<std::endl;
        if (seedSize < 10)
        {
            tAlg.subsimWithTrunc(seedSize, Arg._eps, delta);
        }
        else
        {
            tAlg.subsimWithHIST(seedSize, Arg._eps, delta);
        }
    }
    time(&end);
    cout << "time：" << (end-start) << "秒" << endl; 
    cout << "graph.size()：" << graph.size() << endl; 
    TIO::WriteResult(Arg._outFileName, tRes, Arg._resultFolder);
    TIO::WriteOrderSeeds(Arg._outFileName, tRes, Arg._resultFolder);
    std::cout << "---The End of " << Arg._outFileName << "---\n";
    //analyze the seedsets
    graphLoad->analyzeSeed(Arg._outFileName);
    tAlg.RefreshHypergraph();
    
    return 0;
}