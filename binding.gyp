{
  "targets": [
    {
      "target_name": "towerdef",
      "sources": [  "./cpp/ETP/Point.cpp",
                    "./cpp/index.cpp",
                    "./cpp/clipper/clipper.cpp",
                    "./cpp/core.cpp",
                    "./cpp/utilz.cpp",
                    "./cpp/RVO2/Agent.cpp",
                    "./cpp/RVO2/Definitions.cpp",
                    "./cpp/RVO2/KdTree.cpp",
                    "./cpp/RVO2/Obstacle.cpp",
                    "./cpp/RVO2/RVOSimulator.cpp",
                    "./cpp/RVO2/Vector2.cpp",
                    "./cpp/RVO2/Definitions.cpp",
                    "./cpp/poly2tri/poly2tri.cpp",
                    "./cpp/poly2tri/common/shapes.cc",
                    "./cpp/poly2tri/common/utils.cpp",
                    "./cpp/poly2tri/sweep/advancing_front.cc",
                    "./cpp/poly2tri/sweep/cdt.cc",
                    "./cpp/poly2tri/sweep/sweep_context.cc",
                    "./cpp/poly2tri/sweep/sweep.cc"

      ],
      "include_dirs" : [
          "<!(node -e \"require('nan')\")"
      ],
      'cflags_cc': [ '-fexceptions' ]
    }
  ]
}
