#ifndef _cluster_included_
#define _cluster_included_

bool InVector(vector<int> &v, int index)
{
        int Size=v.size();

        for (int i=0;i<Size;i++)
        {
                if (v[i]==index) return true;
        }
        return false;
}

void AssignCentroidsAtRandom(int Size, int NumClusters, vector<int> &centroids)
{
        int NumPicked=0, pick;
        centroids.clear();
        while (true)
        {
                pick=RandInt(Size);
                if (!InVector(centroids, pick))
                {
                        SafePushBack(centroids, pick, "AssignCentroidsAtRandom");
                        if (centroids.size()>=NumClusters) break;
                }
        }
}

int FindBestCluster(vector<int> &centroids, Matrix &DistanceMatrix, int index)
{
        int NumClusters=centroids.size();
        int BestCluster;
        Real min;

        min=DistanceMatrix[index][centroids[0]];
        BestCluster=0;
        //PrintMatrix(DistanceMatrix);
        for (int i=1;i<NumClusters;i++)
        {
                if (DistanceMatrix[index][centroids[i]]<min)
                {
                        min=DistanceMatrix[index][centroids[i]];
                        BestCluster=i;
                }
        }
        return BestCluster;
}

void AssignClusters(Matrix &DistanceMatrix, vector<int> &centroids, vector< vector<int> > &clusters)
{
        int BestCluster;
        int NumClusters=centroids.size();
        int Size=DistanceMatrix[0].size();
        clusters.clear();
        Safe2DAlloc(clusters, NumClusters, 0, "clusters in AssignClusters");

        for (int i=0;i<Size;i++)
        {
                BestCluster=FindBestCluster(centroids, DistanceMatrix, i);
                SafePushBack(clusters[BestCluster], i, "clusters in AssignClusters");
        }
}

Real CalcSumDistances(int index, vector<int> &cluster, Matrix &DistanceMatrix)
{
        int Size=cluster.size();
        Real SumDistances=0;

        for (int i=0;i<Size;i++)
        {
                SumDistances+=DistanceMatrix[index][cluster[i]];
        }
        return SumDistances;
}

int FindCentroid(Matrix &DistanceMatrix, vector<int> &cluster)
{
        int Size=cluster.size();
        int centroid;
        Real MinSumDistances, SumDistances;

        MinSumDistances=CalcSumDistances(cluster[0], cluster, DistanceMatrix);
        centroid=0;
        for (int i=0;i<Size;i++)
        {
                SumDistances=CalcSumDistances(cluster[i], cluster, DistanceMatrix);
                if (SumDistances<MinSumDistances)
                {
                        MinSumDistances=SumDistances;
                        centroid=cluster[i];
                }
        }
        return centroid;
}

void FindCentroids(Matrix &DistanceMatrix, vector<int> &centroids, vector< vector<int> > &clusters)
{
        int NumClusters=centroids.size();
        int centroid, Size;

        for (int i=0;i<NumClusters;i++)
        {
                centroid=FindCentroid(DistanceMatrix, clusters[i]);
                centroids[i]=centroid;
        }
}

void cluster(Matrix &DistanceMatrix, int NumClusters, vector<int> &centroids, vector< vector<int> > &clusters)
{
        int MaxIterations=1000;
        int Size=DistanceMatrix[0].size();
        AssignCentroidsAtRandom(Size, NumClusters, centroids);
        for (int i=0;i<MaxIterations;i++)
        {
                AssignClusters(DistanceMatrix, centroids, clusters);
                FindCentroids(DistanceMatrix, centroids, clusters);
        }
}

#endif
