{
  "targets": [
    {
      "target_name": "towerdef",
      "sources": [  "./cpp/ETP/Point.cpp",
                    "./cpp/index.cpp",
                    "./cpp/core.cpp",
                    "./cpp/utils.cpp",
                    "./cpp/RVO2/Agent.cpp",
                    "./cpp/RVO2/Definitions.cpp",
                    "./cpp/RVO2/KdTree.cpp",
                    "./cpp/RVO2/Obstacle.cpp",
                    "./cpp/RVO2/RVOSimulator.cpp",
                    "./cpp/RVO2/Vector2.cpp",
                    "./cpp/RVO2/Definitions.cpp",
                    "./cpp/clip2tri-master/clip2tri/clip2tri.cpp",
                    "./cpp/clip2tri-master/clipper/clipper.cpp",
                    "./cpp/clip2tri-master/poly2tri/poly2tri.cpp",
                    "./cpp/clip2tri-master/poly2tri/common/shapes.cc",
                    "./cpp/clip2tri-master/poly2tri/common/utils.cpp",
                    "./cpp/clip2tri-master/poly2tri/sweep/advancing_front.cc",
                    "./cpp/clip2tri-master/poly2tri/sweep/cdt.cc",
                    "./cpp/clip2tri-master/poly2tri/sweep/sweep_context.cc",
                    "./cpp/clip2tri-master/poly2tri/sweep/sweep.cc"

      ],
      "include_dirs" : [
          "<!(node -e \"require('nan')\")"
      ],
      'cflags_cc': [ '-fexceptions' ]
    }
  ]
}
