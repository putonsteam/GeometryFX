//
// Copyright (c) 2016 Advanced Micro Devices, Inc. All rights reserved.
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.
//

#ifndef AMD_GEOMETRYFX_MESH_H
#define AMD_GEOMETRYFX_MESH_H

#include <d3d11.h>
#include <wrl.h>
#include <vector>

#include <DirectXMath.h>

using namespace std;
using namespace DirectX;
namespace AMD
{
namespace GeometryFX_Internal
{

class StaticMesh
{
public:
    StaticMesh(const int vertexCount, const int indexCount, const int meshIndex);

    virtual ~StaticMesh();

public:
    Microsoft::WRL::ComPtr<ID3D11Buffer> vertexBuffer;
    Microsoft::WRL::ComPtr<ID3D11ShaderResourceView> vertexBufferSRV;
    Microsoft::WRL::ComPtr<ID3D11Buffer> indexBuffer;
    Microsoft::WRL::ComPtr<ID3D11ShaderResourceView> indexBufferSRV;
    Microsoft::WRL::ComPtr<ID3D11Buffer> meshConstantsBuffer;

    int vertexCount;
    int faceCount;
    int indexCount;
    int meshIndex;

    int indexOffset;
    int vertexOffset;

    struct Cluster
    {
        DirectX::XMVECTOR aabbMin, aabbMax;
        DirectX::XMVECTOR coneCenter, coneAxis;
        // 面的顺序: 左->右->上->下->前->后
        XMVECTOR frustumPlaneLeft[6];
        XMVECTOR frustumPlaneRight[6];
        XMVECTOR frustumPlaneTop[6];
        XMVECTOR frustumPlaneBottom[6];
        XMVECTOR frustumPlaneFront[6];
        XMVECTOR frustumPlaneBack[6];
        vector<int> indicesLeft;
        vector<int> indicesRight;
        vector<int> indicesTop;
        vector<int> indicesBottom;
        vector<int> indicesFront;
        vector<int> indicesBack;
        vector<int> indices;

        float coneAngleCosine;
        bool valid;
    };

    std::vector<Cluster> clusters;

private:
    StaticMesh(const StaticMesh &);
    StaticMesh &operator=(const StaticMesh &);
};

} // namespace GeometryFX_Internal
} // namespace AMD

#endif // AMD_GEOMETRYFX_MESH_H
