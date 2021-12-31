using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using System;
using System.IO;

public class FVM : MonoBehaviour
{
    float dt = 0.003f;
    float mass = 1;
    float stiffness_0 = 20000.0f;
    float stiffness_1 = 5000.0f;
    float damp = 0.999f;

    int[] Tet;
    int tet_number;         //The number of tetrahedra

    Vector3[] Force;
    Vector3[] V;
    Vector3[] X;
    int number;             //The number of vertices

    Matrix4x4[] inv_Dm;

    //For Laplacian smoothing.
    Vector3[] V_sum;
    int[] V_num;

    SVD svd = new SVD();
    Vector3 floorPosition;
    Vector3 floorNormal;

    const float MuT = 0.5f;
    const float MuN = 0.5f;

    // Start is called before the first frame update
    void Start()
    {
        // FILO IO: Read the house model from files.
        // The model is from Jonathan Schewchuk's Stellar lib.
        {
            string fileContent = File.ReadAllText("Assets/house2.ele");
            string[] Strings = fileContent.Split(new char[] { ' ', '\t', '\r', '\n' }, StringSplitOptions.RemoveEmptyEntries);

            tet_number = int.Parse(Strings[0]);
            Tet = new int[tet_number * 4];

            for (int tet = 0; tet < tet_number; tet++)
            {
                Tet[tet * 4 + 0] = int.Parse(Strings[tet * 5 + 4]) - 1;
                Tet[tet * 4 + 1] = int.Parse(Strings[tet * 5 + 5]) - 1;
                Tet[tet * 4 + 2] = int.Parse(Strings[tet * 5 + 6]) - 1;
                Tet[tet * 4 + 3] = int.Parse(Strings[tet * 5 + 7]) - 1;
            }
        }
        {
            string fileContent = File.ReadAllText("Assets/house2.node");
            string[] Strings = fileContent.Split(new char[] { ' ', '\t', '\r', '\n' }, StringSplitOptions.RemoveEmptyEntries);
            number = int.Parse(Strings[0]);
            X = new Vector3[number];
            for (int i = 0; i < number; i++)
            {
                X[i].x = float.Parse(Strings[i * 5 + 5]) * 0.4f;
                X[i].y = float.Parse(Strings[i * 5 + 6]) * 0.4f;
                X[i].z = float.Parse(Strings[i * 5 + 7]) * 0.4f;
            }
            //Centralize the model.
            Vector3 center = Vector3.zero;
            for (int i = 0; i < number; i++) center += X[i];
            center = center / number;
            for (int i = 0; i < number; i++)
            {
                X[i] -= center;
                float temp = X[i].y;
                X[i].y = X[i].z;
                X[i].z = temp;
            }
        }
        /*tet_number=1;
        Tet = new int[tet_number*4];
        Tet[0]=0;
        Tet[1]=1;
        Tet[2]=2;
        Tet[3]=3;

        number=4;
        X = new Vector3[number];
        V = new Vector3[number];
        Force = new Vector3[number];
        X[0]= new Vector3(0, 0, 0);
        X[1]= new Vector3(1, 0, 0);
        X[2]= new Vector3(0, 1, 0);
        X[3]= new Vector3(0, 0, 1);*/


        //Create triangle mesh.
        Vector3[] vertices = new Vector3[tet_number * 12];
        int vertex_number = 0;
        for (int tet = 0; tet < tet_number; tet++)
        {
            vertices[vertex_number++] = X[Tet[tet * 4 + 0]];
            vertices[vertex_number++] = X[Tet[tet * 4 + 2]];
            vertices[vertex_number++] = X[Tet[tet * 4 + 1]];

            vertices[vertex_number++] = X[Tet[tet * 4 + 0]];
            vertices[vertex_number++] = X[Tet[tet * 4 + 3]];
            vertices[vertex_number++] = X[Tet[tet * 4 + 2]];

            vertices[vertex_number++] = X[Tet[tet * 4 + 0]];
            vertices[vertex_number++] = X[Tet[tet * 4 + 1]];
            vertices[vertex_number++] = X[Tet[tet * 4 + 3]];

            vertices[vertex_number++] = X[Tet[tet * 4 + 1]];
            vertices[vertex_number++] = X[Tet[tet * 4 + 2]];
            vertices[vertex_number++] = X[Tet[tet * 4 + 3]];
        }

        int[] triangles = new int[tet_number * 12];
        for (int t = 0; t < tet_number * 4; t++)
        {
            triangles[t * 3 + 0] = t * 3 + 0;
            triangles[t * 3 + 1] = t * 3 + 1;
            triangles[t * 3 + 2] = t * 3 + 2;
        }
        Mesh mesh = GetComponent<MeshFilter>().mesh;
        mesh.vertices = vertices;
        mesh.triangles = triangles;
        mesh.RecalculateNormals();


        V = new Vector3[number];
        Force = new Vector3[number];
        V_sum = new Vector3[number];
        V_num = new int[number];

        //TODO: Need to allocate and assign inv_Dm

        inv_Dm = new Matrix4x4[tet_number];
        for (int tet = 0; tet < tet_number; tet++)
        {
            inv_Dm[tet] = Build_Edge_Matrix(tet).inverse;
        }

        var floor = GameObject.Find("Floor");
        floorPosition = floor.transform.position;
        floorNormal = new Vector3(0, 1, 0);
        V_sum = new Vector3[number];
        V_num = new int[number];
    }

    Matrix4x4 Build_Edge_Matrix(int tet)
    {
        Matrix4x4 ret = Matrix4x4.zero;
        //TODO: Need to build edge matrix here.
        var x0 = X[Tet[tet * 4]];
        var x1 = X[Tet[tet * 4 + 1]];
        var x2 = X[Tet[tet * 4 + 2]];
        var x3 = X[Tet[tet * 4 + 3]];
        ret.SetColumn(0, x1 - x0);
        ret.SetColumn(1, x2 - x0);
        ret.SetColumn(2, x3 - x0);
        ret.SetColumn(3, new Vector4(0, 0, 0, 1));

        return ret;
    }
    Matrix4x4 MatrixMutipleFloat(Matrix4x4 m, float f)
    {
        Matrix4x4 result = Matrix4x4.zero;
        for (int i = 0; i < 4; i++)
        {
            for (int j = 0; j < 4; j++)
            {
                result[i, j] = m[i, j] * f;
            }
        }
        return result;
    }
    void _Update()
    {
        // Jump up.
        if (Input.GetKeyDown(KeyCode.Space))
        {
            for (int i = 0; i < number; i++)
                V[i].y += 0.2f;
        }

        Vector3 graivity = new Vector3(0, -1, 0) * 9.8f * mass;
        for (int i = 0; i < number; i++)
        {
            //TODO: Add gravity to Force.
            Force[i] = graivity;
            V[i] *= damp;
        }

        for (int tet = 0; tet < tet_number; tet++)
        {
            //TODO: Deformation Gradient
            var F = Build_Edge_Matrix(tet) * inv_Dm[tet];
            //TODO: Green Strain
            var G = F.transpose * F;
            for (int i = 0; i < 4; i++)
            {
                G[i, i] -= 1;
            }
            G = MatrixMutipleFloat(G, 0.5f);
            //TODO: Second PK Stress
            var S = MatrixMutipleFloat(G, 2 * stiffness_1);
            var traceG = G[0, 0] + G[1, 1] + G[2, 2];
            var S2 = MatrixMutipleFloat(Matrix4x4.identity, traceG * stiffness_0);
            for (int i = 0; i < 4; i++)
                for (int j = 0; j < 4; j++)
                    S[i, j] += S2[i, j];
            //TODO: Elastic Force
            var P = F * S;
            var fmatrix = MatrixMutipleFloat(P * inv_Dm[tet].transpose, -1f / 6f / inv_Dm[tet].determinant);
            Vector3[] f = new Vector3[4];
            for (int i = 0; i < 3; i++)
            {
                f[i + 1] = fmatrix.GetColumn(i);
                f[0] -= f[i + 1];
            }
            for (int i = 0; i < 4; i++)
            {
                Force[Tet[tet * 4 + i]] += f[i];
            }
        }
        // laplace
        for (int i = 0; i < number; i ++)
        {
            V[i] += Force[i] / mass * dt;
            V_sum[i] = Vector3.zero;
            V_num[i] = 0;
        }
        for(int tet = 0; tet < tet_number; tet++)
        {
            Vector3 sum = Vector3.zero;
            for(int i = 0; i < 4; i ++)
            {
                sum += V[Tet[tet * 4 + i]];
            }
            for(int i = 0; i < 4; i ++)
            {
                var v_index = Tet[tet * 4 + i];
                V_sum[v_index] += sum - V[v_index];
                V_num[v_index] += 3;
            }
        }
        const float w = 0.5f;
        for (int i = 0; i < number; i++)
        {
            //TODO: Update X and V here.
            V[i] = V[i] * w + V_sum[i] / V_num[i] * (1 - w);
            // I found it's OK just to update x once after collision
            //X[i] += V[i] * dt;
            //TODO: (Particle) collision with floor.
            var phi = X[i] - floorPosition;
            if (Vector3.Dot(phi, floorNormal) < 0)
            {
                if (Vector3.Dot(V[i], floorNormal) < 0)
                {
                    var vNi = Vector3.Dot(V[i], floorNormal) * floorNormal;
                    var vTi = V[i] - vNi;
                    var alpha = Mathf.Max(1 - MuT * (1 + MuT) * vNi.sqrMagnitude / vTi.sqrMagnitude, 0);
                    var vNi_new = -MuN * vNi;
                    var vTi_new = alpha * vTi;
                    V[i] = vNi_new + vTi_new;
                }
            }
            X[i] += V[i] * dt;
        }
    }

    // Update is called once per frame
    void Update()
    {
        for (int l = 0; l < 10; l++)
            _Update();

        // Dump the vertex array for rendering.
        Vector3[] vertices = new Vector3[tet_number * 12];
        int vertex_number = 0;
        for (int tet = 0; tet < tet_number; tet++)
        {
            vertices[vertex_number++] = X[Tet[tet * 4 + 0]];
            vertices[vertex_number++] = X[Tet[tet * 4 + 2]];
            vertices[vertex_number++] = X[Tet[tet * 4 + 1]];
            vertices[vertex_number++] = X[Tet[tet * 4 + 0]];
            vertices[vertex_number++] = X[Tet[tet * 4 + 3]];
            vertices[vertex_number++] = X[Tet[tet * 4 + 2]];
            vertices[vertex_number++] = X[Tet[tet * 4 + 0]];
            vertices[vertex_number++] = X[Tet[tet * 4 + 1]];
            vertices[vertex_number++] = X[Tet[tet * 4 + 3]];
            vertices[vertex_number++] = X[Tet[tet * 4 + 1]];
            vertices[vertex_number++] = X[Tet[tet * 4 + 2]];
            vertices[vertex_number++] = X[Tet[tet * 4 + 3]];
        }
        Mesh mesh = GetComponent<MeshFilter>().mesh;
        mesh.vertices = vertices;
        mesh.RecalculateNormals();
    }
}
