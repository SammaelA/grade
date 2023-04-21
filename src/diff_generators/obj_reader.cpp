#pragma once
#include <vector>
#include "vectors.h"
#include <string>

namespace dgen
{
  std::vector<std::string> explode(std::string aStr, char aDelim)
  {
    std::vector<std::string> res;
    std::string str = aStr.substr(0, aStr.find(aDelim));

    while (str.size() < aStr.size())
    {
      res.push_back(str);
      aStr = aStr.substr(aStr.find(aDelim) + 1);
      str = aStr.substr(0, aStr.find(aDelim));
    }

    res.push_back(str);

    return res;
  }

  std::vector<float> obj2Vec(const std::string aInFilename, bool aVerbose)
  {
    if (aVerbose) { std::cout << "Loading OBJ file <" << aInFilename << ">" << std::endl; }

    // Open file
    std::ifstream objFile(aInFilename.c_str());

    if (objFile.fail())
    {
      std::cout << "Error: could not open file <" << aInFilename << ">" << std::endl;
      exit(1);
    }


    // Extract verts, normals, textures, and faces
    std::vector<float> verts, norms, textures;
    std::vector<float> faces;

    std::vector<float> finalVerts;

    std::string line;

    if (aVerbose) std::cout << "Extracting values from file" << std::endl;

    // Visit each line of the obj file
    while (std::getline(objFile, line))
    {
      // Extract vertex
      // Line starts with v[space]...
      if (line[0] == 'v' && line[1] == ' ')
      {
        std::string lineVals = line.substr(2);
        float val;

        std::string val0 = lineVals.substr(0, lineVals.find(' '));
        val = (float)atof(val0.c_str());
        verts.push_back(val);

        std::string val1 = lineVals.substr(val0.length() + 1,
          lineVals.find(' '));
        val = (float)atof(val1.c_str());
        verts.push_back(val);

        std::string val2 = lineVals.substr(lineVals.find_last_of(' ') + 1);
        val = (float)atof(val2.c_str());
        verts.push_back(val);
      }


      // Extract textures
      // Line starts with vt[space]...
      else if (line[0] == 'v' && line[1] == 't' && line[2] == ' ')
      {
        std::string lineVals = line.substr(3);
        float val;

        std::string val0 = lineVals.substr(0, lineVals.find(' '));
        val = (float)atof(val0.c_str());
        textures.push_back(val);

        std::string val1 = lineVals.substr(val0.length() + 1,
          lineVals.find(' '));
        val = (float)atof(val1.c_str());
        textures.push_back(val);
      }


      // Extract normals
      // Line starts with vn[space]...
      else if (line[0] == 'v' && line[1] == 'n' && line[2] == ' ')
      {
        std::string lineVals = line.substr(3);
        float val;

        std::string val0 = lineVals.substr(0, lineVals.find(' '));
        val = (float)atof(val0.c_str());
        norms.push_back(val);

        std::string val1 = lineVals.substr(val0.length() + 1,
          lineVals.find(' '));
        val = (float)atof(val1.c_str());
        norms.push_back(val);

        std::string val2 = lineVals.substr(lineVals.find_last_of(' ') + 1);
        val = (float)atof(val2.c_str());
        norms.push_back(val);
      }


      //
      // 2. Hash faces
      //
        // Extract faces
        // Line starts with f[space]...
      else if (line[0] == 'f' && line[1] == ' ')
      {
        std::string lineVals = line.substr(2);

        std::string val0 = lineVals.substr(0, lineVals.find_first_of(' '));

        // If the value for this face includes texture and/or 
        // normal, parse them out
        if (val0.find('/') == std::string::npos)
        {
          std::string lineVals = line.substr(2);
          int val;

          std::string val0 = lineVals.substr(0, lineVals.find(' '));
          val = (int)atoi(val0.c_str());
          finalVerts.push_back(verts[(val - 1) * 3]);
          finalVerts.push_back(verts[(val - 1) * 3 + 1]);
          finalVerts.push_back(verts[(val - 1) * 3 + 2]);
          finalVerts.push_back(0);
          finalVerts.push_back(0);
          finalVerts.push_back(0);
          finalVerts.push_back(0);
          finalVerts.push_back(0);

          std::string val1 = lineVals.substr(val0.length() + 1,
            lineVals.find(' '));
          val = (int)atoi(val1.c_str());
          finalVerts.push_back(verts[(val - 1) * 3]);
          finalVerts.push_back(verts[(val - 1) * 3 + 1]);
          finalVerts.push_back(verts[(val - 1) * 3 + 2]);
          finalVerts.push_back(0);
          finalVerts.push_back(0);
          finalVerts.push_back(0);
          finalVerts.push_back(0);
          finalVerts.push_back(0);

          std::string val2 = lineVals.substr(lineVals.find_last_of(' ') + 1);
          val = (int)atoi(val2.c_str());
          finalVerts.push_back(verts[(val - 1) * 3]);
          finalVerts.push_back(verts[(val - 1) * 3 + 1]);
          finalVerts.push_back(verts[(val - 1) * 3 + 2]);
          finalVerts.push_back(0);
          finalVerts.push_back(0);
          finalVerts.push_back(0);
          finalVerts.push_back(0);
          finalVerts.push_back(0);
        }
        else
        {
          // Get first group of values
          std::string g[3];
          g[0] = val0.substr(0, val0.find(' '));

          // Get second group of values
          g[1] = lineVals.substr(line.find(' ') + 1);
          g[1] = g[1].substr(g[1].find(' ') + 1);
          g[1] = g[1].substr(0, g[1].find(' '));

          g[2] = line.substr(line.find_last_of(' ') + 1);

          if (aVerbose)
            std::cout << "Face: (" << g[0] << ") (" << g[1] << ") (" << g[2] << ")" << std::endl;

          for (int k = 0; k < 3; ++k)
          {
            std::vector<std::string> data = explode(g[k], '/');
            if (data[0] != "")
            {
              int val = (int)atoi(data[0].c_str());
              finalVerts.push_back(verts[(val - 1) * 3]);
              finalVerts.push_back(verts[(val - 1) * 3 + 1]);
              finalVerts.push_back(verts[(val - 1) * 3 + 2]);
            }
            else
            {
              finalVerts.push_back(0);
              finalVerts.push_back(0);
              finalVerts.push_back(0);
            }
            if (data[1] != "")
            {
              int val = (int)atoi(data[1].c_str());
              finalVerts.push_back(verts[(val - 1) * 3]);
              finalVerts.push_back(verts[(val - 1) * 3 + 1]);
              finalVerts.push_back(verts[(val - 1) * 3 + 2]);
            }
            else
            {
              finalVerts.push_back(0);
              finalVerts.push_back(0);
              finalVerts.push_back(0);
            }
            if (data.size() > 2 && data[2] != "")
            {
              int val = (int)atoi(data[2].c_str());
              finalVerts.push_back(verts[(val - 1) * 2]);
              finalVerts.push_back(verts[(val - 1) * 2 + 1]);
            }
            else
            {
              finalVerts.push_back(0);
              finalVerts.push_back(0);
            }
          }

        }
      }
    } /* end getline(file, line) */

    if (aVerbose) std::cout << "Finished extracting values from file" << std::endl
      << "Quick count check:" << std::endl
      << "\tVerts = " << verts.size() << std::endl
      << "\tNorms = " << norms.size() << std::endl
      << "\tTexts = " << textures.size() << std::endl
      << "\tFaces = " << finalVerts.size() << std::endl;

    objFile.close();

    if (aVerbose) std::cout << "Preparing to build faces" << std::endl;

    return finalVerts;
    
  }
}