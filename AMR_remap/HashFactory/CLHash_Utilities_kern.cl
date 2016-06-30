static inline unsigned int intintHash_CompressIdentity(char data, int hashCode){ return hashCode; }
typedef struct intintHash_CompressLCGData{ 
    long unsigned int a;
    long unsigned int c;
    unsigned int m;
    unsigned int n;
}intintHash_CompressLCGData;

static inline unsigned int intintHash_CompressLCG(intintHash_CompressLCGData compressLCGData, int hashCode){
    return ((compressLCGData.a * hashCode + compressLCGData.c) % compressLCGData.m) % compressLCGData.n; 
}

typedef struct intintLCGQuadraticOpenCompactCLHash_TableData{ 
    int hashID;
    unsigned int numBuckets;
    intintHash_CompressLCGData compressFuncData;
}intintLCGQuadraticOpenCompactCLHash_TableData; 

typedef struct intintLCGQuadraticOpenCompactCLHash_Bucket{
    int key; 
    int value; 
}intintLCGQuadraticOpenCompactCLHash_Bucket; 

__kernel void intintLCGQuadraticOpenCompactCLHash_Empty(__global char *tableData){ 
    uint index = get_global_id(0); 
    if(index >= ((__global intintLCGQuadraticOpenCompactCLHash_TableData*)tableData)->numBuckets){ 
        return; 
    } 
    __global intintLCGQuadraticOpenCompactCLHash_Bucket *buckets 
    = (__global intintLCGQuadraticOpenCompactCLHash_Bucket*)
      &tableData[sizeof(intintLCGQuadraticOpenCompactCLHash_TableData)]; 
      buckets[index].key = UINT_MAX;
}

typedef struct intintIdentityPerfectCLHash_TableData{ 
    int hashID; 
    unsigned int numBuckets; 
    char compressFuncData;
}intintIdentityPerfectCLHash_TableData; 

typedef struct intintIdentityPerfectCLHash_Bucket{ 
    int key; 
    int value; 
}intintIdentityPerfectCLHash_Bucket; 

__kernel void intintIdentityPerfectCLHash_Empty(__global char *tableData){ 
    uint index = get_global_id(0); 
    if(index >= ((__global intintIdentityPerfectCLHash_TableData*)tableData)->numBuckets){ return; } 
    __global intintIdentityPerfectCLHash_Bucket *buckets 
      = (__global intintIdentityPerfectCLHash_Bucket*)
        &tableData[sizeof(intintIdentityPerfectCLHash_TableData)]; 
    buckets[index].key = UINT_MAX;
}

typedef struct intintIdentitySentinelPerfectCLHash_TableData{ 
    int hashID; 
    unsigned int numBuckets; 
    char compressFuncData; 
    uint emptyValue; 
}intintIdentitySentinelPerfectCLHash_TableData;
 
typedef struct intintIdentitySentinelPerfectCLHash_Bucket{ uint value; }
intintIdentitySentinelPerfectCLHash_Bucket; 

__kernel void intintIdentitySentinelPerfectCLHash_Empty(__global char *tableData){ 
    uint index = get_global_id(0); 
    if(index >= ((__global intintIdentitySentinelPerfectCLHash_TableData*)tableData)->numBuckets){ 
        return;
    }
    __global intintIdentitySentinelPerfectCLHash_Bucket *buckets 
      = (__global intintIdentitySentinelPerfectCLHash_Bucket*)
        &tableData[sizeof(intintIdentitySentinelPerfectCLHash_TableData)];
    buckets[index].value 
      = ((__global intintIdentitySentinelPerfectCLHash_TableData*)tableData)->emptyValue; 
}

typedef struct intintLCGLinearOpenCompactCLHash_TableData{ 
    int hashID; 
    unsigned int numBuckets; 
    intintHash_CompressLCGData compressFuncData;
}intintLCGLinearOpenCompactCLHash_TableData; 

typedef struct intintLCGLinearOpenCompactCLHash_Bucket{ 
    int key; 
    int value; 
}intintLCGLinearOpenCompactCLHash_Bucket; 

__kernel void intintLCGLinearOpenCompactCLHash_Empty(__global char *tableData){ 
    uint index = get_global_id(0); 
    if(index >= ((__global intintLCGLinearOpenCompactCLHash_TableData*)tableData)->numBuckets){ 
        return; 
    } 
    __global intintLCGLinearOpenCompactCLHash_Bucket *buckets 
      = (__global intintLCGLinearOpenCompactCLHash_Bucket*)
        &tableData[sizeof(intintLCGLinearOpenCompactCLHash_TableData)];
    buckets[index].key = UINT_MAX;
}
