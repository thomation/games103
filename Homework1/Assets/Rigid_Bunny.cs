﻿using UnityEngine;
using System.Collections;

public class Rigid_Bunny : MonoBehaviour 
{
	bool launched 		= false;
	float dt 			= 0.015f;
	Vector3 v 			= new Vector3(0, 0, 0);	// velocity
	Vector3 w 			= new Vector3(0, 0, 0);	// angular velocity
	
	float mass;									// mass
	Matrix4x4 I_ref;							// reference inertia

	float linear_decay	= 0.999f;				// for velocity decay
	float angular_decay	= 0.98f;				
	float MuT = 0.5f;
	float MuN = 0.5f;


	// Use this for initialization
	void Start () 
	{		
		Mesh mesh = GetComponent<MeshFilter>().mesh;
		Vector3[] vertices = mesh.vertices;

		float m=1;
		mass=0;
		for (int i=0; i<vertices.Length; i++) 
		{
			mass += m;
			float diag=m*vertices[i].sqrMagnitude;
			I_ref[0, 0]+=diag;
			I_ref[1, 1]+=diag;
			I_ref[2, 2]+=diag;
			I_ref[0, 0]-=m*vertices[i][0]*vertices[i][0];
			I_ref[0, 1]-=m*vertices[i][0]*vertices[i][1];
			I_ref[0, 2]-=m*vertices[i][0]*vertices[i][2];
			I_ref[1, 0]-=m*vertices[i][1]*vertices[i][0];
			I_ref[1, 1]-=m*vertices[i][1]*vertices[i][1];
			I_ref[1, 2]-=m*vertices[i][1]*vertices[i][2];
			I_ref[2, 0]-=m*vertices[i][2]*vertices[i][0];
			I_ref[2, 1]-=m*vertices[i][2]*vertices[i][1];
			I_ref[2, 2]-=m*vertices[i][2]*vertices[i][2];
		}
		I_ref [3, 3] = 1;
	}
    Matrix4x4 Get_Cross_Matrix(Vector3 a)
	{
		//Get the cross product matrix of vector a
		Matrix4x4 A = Matrix4x4.zero;
		A [0, 0] = 0; 
		A [0, 1] = -a [2]; 
		A [0, 2] = a [1]; 
		A [1, 0] = a [2]; 
		A [1, 1] = 0; 
		A [1, 2] = -a [0]; 
		A [2, 0] = -a [1]; 
		A [2, 1] = a [0]; 
		A [2, 2] = 0; 
		A [3, 3] = 1;
		return A;
	}

	// In this function, update v and w by the impulse due to the collision with
	//a plane <P, N>
	void Collision_Impulse(Vector3 P, Vector3 N)
	{
		var insidePoint = Vector3.zero;
		var insideAmount = 0;
		Mesh mesh = GetComponent<MeshFilter>().mesh;
		Vector3[] vertices = mesh.vertices;
		Matrix4x4 R = Matrix4x4.Rotate(transform.rotation);
        for (int i = 0; i < vertices.Length; i++)
        {
			Vector3 Rri = R.MultiplyVector(vertices[i]);
            Vector3 xi = transform.position + Rri;
            var d = Vector3.Dot(xi - P, N);
			if (d < 0)
			{
				insideAmount++;
				insidePoint += vertices[i];
			}
        }
        if (insideAmount > 0)
        {
            insidePoint /= insideAmount;
			Vector3 Rri = R.MultiplyVector(insidePoint);
            var vi = v + Vector3.Cross(w, Rri);
            if (Vector3.Dot(vi, N) < 0)
            {
				var vNi = Vector3.Dot(vi, N) * N;
				var vTi = vi - vNi;
				var alpha = Mathf.Max(1 - MuT * (1 + MuT) * vNi.sqrMagnitude / vTi.sqrMagnitude, 0);
				var vNi_new = -MuN * vNi;
				var vTi_new = alpha * vTi;
				var vi_new = vNi_new + vTi_new;

				var K = Matrix4x4.zero;
				K[0, 0] = K[1, 1] = K[2, 2] = 1f / mass;
				var Rri_cross = Get_Cross_Matrix(Rri);
				var I_ref_inverse = I_ref.inverse;
				var sub = Rri_cross * I_ref_inverse * Rri_cross;
				for(int ri = 0; ri < 4; ri ++)
					for(int ci = 0; ci < 4; ci ++)
                        K[ri, ci] = K[ri, ci] - sub[ri, ci];

				var j = K.inverse.MultiplyVector(vi_new - vi);
				v += j / mass;
				w += I_ref_inverse.MultiplyVector(Vector3.Cross(Rri, j));
            }
        }
    }

    // Update is called once per frame
    void Update () 
	{
		//Game Control
		if(Input.GetKey("r"))
		{
			transform.position = new Vector3 (0, 0.6f, 0);
            transform.rotation = new Quaternion(0, 0, 0, 0);
			launched=false;
		}
		if(Input.GetKey("l"))
		{
			v = new Vector3 (5, 2, 0);
			w = new Vector3(0, 50, 0);
			launched=true;
		}
		if (!launched)
			return;
		// Part I: Update velocities
		var f = mass * new Vector3(0, -9.8f, 0);
		v += f * dt / mass;
		v *= linear_decay;

		w *= angular_decay;


        // Part II: Collision Impulse
        Collision_Impulse(new Vector3(0, 0.01f, 0), new Vector3(0, 1, 0));
		Collision_Impulse(new Vector3(2, 0, 0), new Vector3(-1, 0, 0));

		// Part III: Update position & orientation
		//Update linear status
		Vector3 x    = transform.position;
		x += v * dt;
		//Update angular status
		Quaternion q = transform.rotation;
		var dw = w * dt * 0.5f;
		var dq = new Quaternion(dw.x, dw.y, dw.z, 0f) * q;
		q = new Quaternion(q.x + dq.x, q.y + dq.y, q.z + dq.z, q.w + dq.w);
		// Part IV: Assign to the object
		transform.position = x;
		transform.rotation = q;
	}
}
