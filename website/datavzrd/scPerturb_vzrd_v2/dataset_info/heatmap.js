const heatmap_plot = {
    "$schema": "https://vega.github.io/schema/vega-lite/v5.json",
    "data": {
        "values": []
    },
    "hconcat": [
        
        {
            
            "layer": [
            
            {
                "mark": {
                    "type": "text",
                    "align":"left"
                },
                "encoding": {
                    "x": {
                        "field": "dummy",
                        "axis": {
                            "labelExpr": "'Download Link'",
                            "labelAngle": -90,
                            "labelOffset": 7,
                            "title": null,
                            "ticks": false,
                            "orient": "top",
                            "domain": false
                        }
                    },
                    "y": {
                        "field": "index",
                        "type": "ordinal",
                        "axis": null,
                    },
                    "text": {
                        "field": "Download Link",
                        "type": "nominal",
                        "title": "Download Link",
                        
                    }
                }
            }
            
        ]},
        
        {
            
            "layer": [
            
            {
                "mark": {
                    "type": "text",
                    "align":"left"
                },
                "encoding": {
                    "x": {
                        "field": "dummy",
                        "axis": {
                            "labelExpr": "'Publication index'",
                            "labelAngle": -90,
                            "labelOffset": 7,
                            "title": null,
                            "ticks": false,
                            "orient": "top",
                            "domain": false
                        }
                    },
                    "y": {
                        "field": "index",
                        "type": "ordinal",
                        "axis": null,
                    },
                    "text": {
                        "field": "Publication index",
                        "type": "nominal",
                        "title": "Publication index",
                        
                    }
                }
            }
            
        ]},
        
        {
            
            "layer": [
            
            {
                "mark": {
                    "type": "text",
                    "align":"left"
                },
                "encoding": {
                    "x": {
                        "field": "dummy",
                        "axis": {
                            "labelExpr": "'Dataset index'",
                            "labelAngle": -90,
                            "labelOffset": 7,
                            "title": null,
                            "ticks": false,
                            "orient": "top",
                            "domain": false
                        }
                    },
                    "y": {
                        "field": "index",
                        "type": "ordinal",
                        "axis": null,
                    },
                    "text": {
                        "field": "Dataset index",
                        "type": "nominal",
                        "title": "Dataset index",
                        
                    }
                }
            }
            
        ]},
        
        {
            
            "layer": [
            
            {
                "mark": {
                    "type": "text",
                    "align":"left"
                },
                "encoding": {
                    "x": {
                        "field": "dummy",
                        "axis": {
                            "labelExpr": "'Title'",
                            "labelAngle": -90,
                            "labelOffset": 7,
                            "title": null,
                            "ticks": false,
                            "orient": "top",
                            "domain": false
                        }
                    },
                    "y": {
                        "field": "index",
                        "type": "ordinal",
                        "axis": null,
                    },
                    "text": {
                        "field": "Title",
                        "type": "nominal",
                        "title": "Title",
                        
                    }
                }
            }
            
        ]},
        
        {
            
            "layer": [
            
            {
                "mark": {
                    "type": "text",
                    "align":"left"
                },
                "encoding": {
                    "x": {
                        "field": "dummy",
                        "axis": {
                            "labelExpr": "'First Author'",
                            "labelAngle": -90,
                            "labelOffset": 7,
                            "title": null,
                            "ticks": false,
                            "orient": "top",
                            "domain": false
                        }
                    },
                    "y": {
                        "field": "index",
                        "type": "ordinal",
                        "axis": null,
                    },
                    "text": {
                        "field": "First Author",
                        "type": "nominal",
                        "title": "First Author",
                        
                    }
                }
            }
            
        ]},
        
        {
            
            "layer": [
            
            {
                "mark": {
                    "type": "text",
                    "align":"left"
                },
                "encoding": {
                    "x": {
                        "field": "dummy",
                        "axis": {
                            "labelExpr": "'Organisms'",
                            "labelAngle": -90,
                            "labelOffset": 7,
                            "title": null,
                            "ticks": false,
                            "orient": "top",
                            "domain": false
                        }
                    },
                    "y": {
                        "field": "index",
                        "type": "ordinal",
                        "axis": null,
                    },
                    "text": {
                        "field": "Organisms",
                        "type": "nominal",
                        "title": "Organisms",
                        
                    }
                }
            }
            
        ]},
        
        {
            
            "layer": [
            
            {
                "mark": {
                    "type": "rect",
                    "align":"left",
                    "width": 138
                },
                "encoding": {
                    "x": {
                        "field": "dummy",
                        "axis": {
                            "labelExpr": "'Modality'",
                            "labelAngle": -90,
                            
                            "title": null,
                            "ticks": false,
                            "orient": "top",
                            "domain": false
                        }
                    },
                    "y": {
                        "field": "index",
                        "type": "ordinal",
                        "axis": null,
                    },
                    "color": {
                        "field": "Modality",
                        "type": "nominal",
                        "title": "Modality",
                        "scale": {
                            "domain": ["RNA","RNA + protein (RNA)","ATAC","ATAC + protein","RNA + protein (protein)"],
                            "type": "ordinal",
                            
                            "range": ["#1f77b4","#1f77b4","#ff7f0e","#bcbd22","#d62728"],
                            
                        },
                        "legend": null
                    }
                }
            }
            ,
            {
                "mark": {
                    "type": "text",
                    "align":"left",
                    "xOffset": -6
                },
                "encoding": {
                    "text": {"field": "Modality"},
                    "y": {
                        "field": "index",
                        "type": "ordinal",
                        "axis": null,
                    }
                }
            }
        
        ]},
        
        {
            
            "layer": [
            
            {
                "mark": {
                    "type": "text",
                    "align":"left"
                },
                "encoding": {
                    "x": {
                        "field": "dummy",
                        "axis": {
                            "labelExpr": "'Method'",
                            "labelAngle": -90,
                            "labelOffset": 7,
                            "title": null,
                            "ticks": false,
                            "orient": "top",
                            "domain": false
                        }
                    },
                    "y": {
                        "field": "index",
                        "type": "ordinal",
                        "axis": null,
                    },
                    "text": {
                        "field": "Method",
                        "type": "nominal",
                        "title": "Method",
                        
                    }
                }
            }
            
        ]},
        
        {
            
            "layer": [
            
            {
                "mark": {
                    "type": "text",
                    "align":"left"
                },
                "encoding": {
                    "x": {
                        "field": "dummy",
                        "axis": {
                            "labelExpr": "'Tissues'",
                            "labelAngle": -90,
                            "labelOffset": 7,
                            "title": null,
                            "ticks": false,
                            "orient": "top",
                            "domain": false
                        }
                    },
                    "y": {
                        "field": "index",
                        "type": "ordinal",
                        "axis": null,
                    },
                    "text": {
                        "field": "Tissues",
                        "type": "nominal",
                        "title": "Tissues",
                        
                    }
                }
            }
            
        ]},
        
        {
            
            "layer": [
            
            {
                "mark": {
                    "type": "text",
                    "align":"left"
                },
                "encoding": {
                    "x": {
                        "field": "dummy",
                        "axis": {
                            "labelExpr": "'Perturbation'",
                            "labelAngle": -90,
                            "labelOffset": 7,
                            "title": null,
                            "ticks": false,
                            "orient": "top",
                            "domain": false
                        }
                    },
                    "y": {
                        "field": "index",
                        "type": "ordinal",
                        "axis": null,
                    },
                    "text": {
                        "field": "Perturbation",
                        "type": "nominal",
                        "title": "Perturbation",
                        
                    }
                }
            }
            
        ]},
        
        {
            
            "layer": [
            
            {
                "mark": {
                    "type": "text",
                    "align":"left"
                },
                "encoding": {
                    "x": {
                        "field": "dummy",
                        "axis": {
                            "labelExpr": "'Disease'",
                            "labelAngle": -90,
                            "labelOffset": 7,
                            "title": null,
                            "ticks": false,
                            "orient": "top",
                            "domain": false
                        }
                    },
                    "y": {
                        "field": "index",
                        "type": "ordinal",
                        "axis": null,
                    },
                    "text": {
                        "field": "Disease",
                        "type": "nominal",
                        "title": "Disease",
                        
                    }
                }
            }
            
        ]},
        
        {
            
            "layer": [
            
            {
                "mark": {
                    "type": "text",
                    "align":"left"
                },
                "encoding": {
                    "x": {
                        "field": "dummy",
                        "axis": {
                            "labelExpr": "'Cell Type'",
                            "labelAngle": -90,
                            "labelOffset": 7,
                            "title": null,
                            "ticks": false,
                            "orient": "top",
                            "domain": false
                        }
                    },
                    "y": {
                        "field": "index",
                        "type": "ordinal",
                        "axis": null,
                    },
                    "text": {
                        "field": "Cell Type",
                        "type": "nominal",
                        "title": "Cell Type",
                        
                    }
                }
            }
            
        ]},
        
        {
            
            "layer": [
            
            {
                "mark": {
                    "type": "text",
                    "align":"left"
                },
                "encoding": {
                    "x": {
                        "field": "dummy",
                        "axis": {
                            "labelExpr": "'Tissue Type'",
                            "labelAngle": -90,
                            "labelOffset": 7,
                            "title": null,
                            "ticks": false,
                            "orient": "top",
                            "domain": false
                        }
                    },
                    "y": {
                        "field": "index",
                        "type": "ordinal",
                        "axis": null,
                    },
                    "text": {
                        "field": "Tissue Type",
                        "type": "nominal",
                        "title": "Tissue Type",
                        
                    }
                }
            }
            
        ]},
        
        {
            
            "layer": [
            
            {
                "mark": {
                    "type": "rect",
                    "align":"left",
                    "width": 20
                },
                "encoding": {
                    "x": {
                        "field": "dummy",
                        "axis": {
                            "labelExpr": "'Cancer'",
                            "labelAngle": -90,
                            
                            "title": null,
                            "ticks": false,
                            "orient": "top",
                            "domain": false
                        }
                    },
                    "y": {
                        "field": "index",
                        "type": "ordinal",
                        "axis": null,
                    },
                    "color": {
                        "field": "Cancer",
                        "type": "nominal",
                        "title": "Cancer",
                        "scale": {
                            "domain": ["y","n"],
                            "type": "ordinal",
                            
                            "range": ["#2CA02C","#D62728"],
                            
                        },
                        "legend": null
                    }
                }
            }
            ,
            {
                "mark": {
                    "type": "text",
                    "align":"left",
                    "xOffset": -6
                },
                "encoding": {
                    "text": {"field": "Cancer"},
                    "y": {
                        "field": "index",
                        "type": "ordinal",
                        "axis": null,
                    }
                }
            }
        
        ]},
        
        {
            
            "layer": [
            
            {
                "mark": {
                    "type": "rect",
                    "align":"left",
                    "width": 24
                },
                "encoding": {
                    "x": {
                        "field": "dummy",
                        "axis": {
                            "labelExpr": "'Year'",
                            "labelAngle": -90,
                            
                            "title": null,
                            "ticks": false,
                            "orient": "top",
                            "domain": false
                        }
                    },
                    "y": {
                        "field": "index",
                        "type": "ordinal",
                        "axis": null,
                    },
                    "color": {
                        "field": "Year",
                        "type": "quantitative",
                        "title": "Year",
                        "scale": {
                            "domain": [2016.0,2023.0],
                            "type": "linear",
                            "scheme": "blues",
                            
                            
                        },
                        "legend": null
                    }
                }
            }
            ,
            {
                "mark": {
                    "type": "text",
                    "align":"left",
                    "xOffset": -6
                },
                "encoding": {
                    "text": {"field": "Year"},
                    "y": {
                        "field": "index",
                        "type": "ordinal",
                        "axis": null,
                    }
                }
            }
        
        ]},
        
        {
            
            "layer": [
            
            {
                "mark": {
                    "type": "text",
                    "align":"left"
                },
                "encoding": {
                    "x": {
                        "field": "dummy",
                        "axis": {
                            "labelExpr": "'Total Number of Cells'",
                            "labelAngle": -90,
                            "labelOffset": 7,
                            "title": null,
                            "ticks": false,
                            "orient": "top",
                            "domain": false
                        }
                    },
                    "y": {
                        "field": "index",
                        "type": "ordinal",
                        "axis": null,
                    },
                    "text": {
                        "field": "Total Number of Cells",
                        "type": "quantitative",
                        "title": "Total Number of Cells",
                        
                    }
                }
            }
            
        ]},
        
        {
            
            "layer": [
            
            {
                "mark": {
                    "type": "text",
                    "align":"left"
                },
                "encoding": {
                    "x": {
                        "field": "dummy",
                        "axis": {
                            "labelExpr": "'Total Number of Perturbations'",
                            "labelAngle": -90,
                            "labelOffset": 7,
                            "title": null,
                            "ticks": false,
                            "orient": "top",
                            "domain": false
                        }
                    },
                    "y": {
                        "field": "index",
                        "type": "ordinal",
                        "axis": null,
                    },
                    "text": {
                        "field": "Total Number of Perturbations",
                        "type": "quantitative",
                        "title": "Total Number of Perturbations",
                        
                    }
                }
            }
            
        ]},
        
        {
            
            "layer": [
            
            {
                "mark": {
                    "type": "text",
                    "align":"left"
                },
                "encoding": {
                    "x": {
                        "field": "dummy",
                        "axis": {
                            "labelExpr": "'Median Number of Genes Expressed per Cell'",
                            "labelAngle": -90,
                            "labelOffset": 7,
                            "title": null,
                            "ticks": false,
                            "orient": "top",
                            "domain": false
                        }
                    },
                    "y": {
                        "field": "index",
                        "type": "ordinal",
                        "axis": null,
                    },
                    "text": {
                        "field": "Median Number of Genes Expressed per Cell",
                        "type": "quantitative",
                        "title": "Median Number of Genes Expressed per Cell",
                        
                    }
                }
            }
            
        ]},
        
        {
            
            "layer": [
            
            {
                "mark": {
                    "type": "text",
                    "align":"left"
                },
                "encoding": {
                    "x": {
                        "field": "dummy",
                        "axis": {
                            "labelExpr": "'Median Number of UMI Counts per Cell'",
                            "labelAngle": -90,
                            "labelOffset": 7,
                            "title": null,
                            "ticks": false,
                            "orient": "top",
                            "domain": false
                        }
                    },
                    "y": {
                        "field": "index",
                        "type": "ordinal",
                        "axis": null,
                    },
                    "text": {
                        "field": "Median Number of UMI Counts per Cell",
                        "type": "quantitative",
                        "title": "Median Number of UMI Counts per Cell",
                        
                    }
                }
            }
            
        ]},
        
        {
            
            "layer": [
            
            {
                "mark": {
                    "type": "text",
                    "align":"left"
                },
                "encoding": {
                    "x": {
                        "field": "dummy",
                        "axis": {
                            "labelExpr": "'Highest Number of Perturbations in one Cell'",
                            "labelAngle": -90,
                            "labelOffset": 7,
                            "title": null,
                            "ticks": false,
                            "orient": "top",
                            "domain": false
                        }
                    },
                    "y": {
                        "field": "index",
                        "type": "ordinal",
                        "axis": null,
                    },
                    "text": {
                        "field": "Highest Number of Perturbations in one Cell",
                        "type": "quantitative",
                        "title": "Highest Number of Perturbations in one Cell",
                        
                    }
                }
            }
            
        ]},
        
        {
            
            "layer": [
            
            {
                "mark": {
                    "type": "text",
                    "align":"left"
                },
                "encoding": {
                    "x": {
                        "field": "dummy",
                        "axis": {
                            "labelExpr": "'Median Number of Cells per Perturbation'",
                            "labelAngle": -90,
                            "labelOffset": 7,
                            "title": null,
                            "ticks": false,
                            "orient": "top",
                            "domain": false
                        }
                    },
                    "y": {
                        "field": "index",
                        "type": "ordinal",
                        "axis": null,
                    },
                    "text": {
                        "field": "Median Number of Cells per Perturbation",
                        "type": "quantitative",
                        "title": "Median Number of Cells per Perturbation",
                        
                    }
                }
            }
            
        ]},
        
        {
            
            "layer": [
            
            {
                "mark": {
                    "type": "text",
                    "align":"left"
                },
                "encoding": {
                    "x": {
                        "field": "dummy",
                        "axis": {
                            "labelExpr": "'Smallest Number of Cells per Perturbation'",
                            "labelAngle": -90,
                            "labelOffset": 7,
                            "title": null,
                            "ticks": false,
                            "orient": "top",
                            "domain": false
                        }
                    },
                    "y": {
                        "field": "index",
                        "type": "ordinal",
                        "axis": null,
                    },
                    "text": {
                        "field": "Smallest Number of Cells per Perturbation",
                        "type": "quantitative",
                        "title": "Smallest Number of Cells per Perturbation",
                        
                    }
                }
            }
            
        ]},
        
        {
            
            "layer": [
            
            {
                "mark": {
                    "type": "text",
                    "align":"left"
                },
                "encoding": {
                    "x": {
                        "field": "dummy",
                        "axis": {
                            "labelExpr": "'Highest Number of Cells per Perturbation'",
                            "labelAngle": -90,
                            "labelOffset": 7,
                            "title": null,
                            "ticks": false,
                            "orient": "top",
                            "domain": false
                        }
                    },
                    "y": {
                        "field": "index",
                        "type": "ordinal",
                        "axis": null,
                    },
                    "text": {
                        "field": "Highest Number of Cells per Perturbation",
                        "type": "quantitative",
                        "title": "Highest Number of Cells per Perturbation",
                        
                    }
                }
            }
            
        ]}
        
    ],
    "config": {
        "style": {"cell": {"stroke": "transparent"}, "guide-label": {"fontWeight": "bold"}},
        "concat": {"spacing": 0},
        "text": {"limit": 135}
    }
}