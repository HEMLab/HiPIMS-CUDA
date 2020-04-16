template<typename T_a, typename T_b>
class functor{
    
public:
    
    inline __host__ __device__
    float operator() (const T_a &x, const T_b &y){
        
        float result = x - y;
        
        return result;
        
    }
    
};




template <typename T>
__global__ void kernel(float *a, float *b, float *c, int N, T f){
    
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    if (i < N){
        float a = 1;
        float b = 2;
        c[i] = f(a,b);
    }
    
}

template <typename T >
void foo(float *a, float *b, float *c, int N, T f){
    
    kernel<<<128,128>>>(a, b, c, N, f);
    
}
