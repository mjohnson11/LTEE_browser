// Utility funciton
function format_name(p) {
  if (p=='Anc') return 'Anc';
  return p[0]=='m' ? 'Ara–'+p[1] : 'Ara+'+p[1];
}

// global variables
const genome_len = 4629812;
let focal_pop = 'Ara–1';
const all_pops = ['Anc', 'm1', 'm2', 'm3', 'm4', 'm5', 'm6', 'p1', 'p2', 'p3', 'p4', 'p5', 'p6'];
const formatted_pops = all_pops.map(format_name);

const allele_graphs = [];
const corr_graph_overlays = []

let search_div;
let tooltip;
let gene_descrip_text;

let sbg;
let spagplot;
let corr_data;
const spaghetti_data = []; // for spaghetti plot across populations
let locusID_list;
let gene_list;
let count_data;
let count_data_map; // map from locusID_list index to rows
const allele_data = {};
const allele_data_map = {};
const freq_data = {};
const freq_data_map = {};
const pop_overlays = [];

let row_clicked;

// layout
const layout = {}
layout.basic_h = 220;
layout.small_h = 70;
layout.small_buf = 5;
layout.main_col_w = 700;
layout.right_col_w = 300;
layout.corr_w = 210;
layout.sgb_h = 110;
layout.buf = 20;
layout.inner_buf = 40;
layout.tn_x = 100;
layout.tn_buf = 30;

layout.sgb_left = 0;
layout.sgb_top = 0;
layout.sgb_width = layout.main_col_w+layout.right_col_w;
layout.sgb_height = layout.sgb_h;

layout.corr_div_left = 0;
layout.corr_div_top = layout.sgb_h+layout.basic_h+layout.buf-5;
layout.corr_div_width = layout.main_col_w;
layout.corr_div_height = layout.basic_h*2;

layout.fit_div_left = layout.main_col_w;
layout.fit_div_top = layout.sgb_h;
layout.fit_div_width = layout.right_col_w;
layout.fit_div_height = layout.basic_h;

layout.focal_allele_left = 0; 
layout.focal_allele_top = layout.sgb_h; 
layout.focal_allele_width = layout.main_col_w; 
layout.focal_allele_height = layout.basic_h;

layout.small_allele_left = layout.main_col_w;
layout.small_allele_top = layout.sgb_h+235;
layout.small_allele_width = layout.right_col_w;



//colors
const hover_color = '#FF4444';
const click_color = '#4444FF';
const base_color = '#FCF0FC';
const background_color = '#000A00';

// axes etc.
const gens_domain = [0, 62000];
const corr_domain = [-0.6, 0.15];
const color_domain = [-0.3, 0.05];
const fit_domain = [0, 0.5];
const gens_x = d3.scaleLinear().range([layout.inner_buf, layout.main_col_w-layout.inner_buf]).domain(gens_domain);
const right_gens_x = d3.scaleLinear().range([layout.inner_buf, layout.right_col_w-layout.inner_buf]).domain(gens_domain);
const right_halfgens_x = d3.scaleLinear().range([layout.small_buf, layout.right_col_w/2-layout.small_buf]).domain(gens_domain);
const fit_y = d3.scaleLinear().range([layout.basic_h-layout.inner_buf, layout.inner_buf]).domain(fit_domain);
const freq_y = d3.scaleLinear().range([layout.basic_h-layout.inner_buf, layout.inner_buf]).domain([0, 1]);
const s_traj = d3.line().x(function(d) { return gens_x(d.x); }).y(function(d) { return fit_y(d.y); });
const corr_x1 = d3.scaleLinear().range([layout.inner_buf, layout.corr_w]).domain(corr_domain);
const corr_y = d3.scaleLinear().range([layout.basic_h-layout.inner_buf, layout.inner_buf]).domain(corr_domain);
const corr_x2 = d3.scaleLinear().range([layout.inner_buf+layout.corr_w, layout.corr_w*2]).domain(corr_domain);
const corr_x3 = d3.scaleLinear().range([layout.inner_buf*2+layout.corr_w*2, layout.corr_w*3+layout.inner_buf]).domain(corr_domain);
const pop_y = d3.scaleBand().range([layout.basic_h*2-layout.inner_buf, layout.inner_buf+layout.basic_h]).domain(formatted_pops);
const pop_y_offset = (pop_y.range()[0]-pop_y.range()[1])/(2*all_pops.length);
const pop_y_for_WOW = d3.scaleLinear().range([layout.basic_h-layout.inner_buf*2, 0]).domain([layout.basic_h*2-layout.inner_buf, layout.inner_buf+layout.basic_h]);
const corr_x3_canvas = d3.scaleLinear().range([0, layout.corr_w-layout.inner_buf]).domain(corr_domain);
const pop_y_array = formatted_pops.map((p) => pop_y(p)+pop_y_offset);

const log_freq_y = d3.scaleLinear().range([layout.basic_h*2-layout.inner_buf, layout.basic_h+layout.inner_buf]).domain([-6, 0]);
const tn_x_scales = [];

const color_scale = d3.scaleDiverging().domain([color_domain[0], 0, color_domain[1]]).interpolator(d3.interpolatePuOr);

// UTILITY AXES FUNCTION //

function standard_axes(svg_obj, xscale, yscale, xlabel, ylabel, title='') {
  const x_ax = d3.axisBottom(xscale).ticks(5);
  const y_ax = d3.axisLeft(yscale);
  if (yscale.domain().length == 2) y_ax.ticks(5);
  const xmid = (xscale.range()[0]+xscale.range()[1])/2;

  const y_called = svg_obj.append('g')
    .attr('class', 'axis')
    .attr('transform', 'translate(' + xscale.range()[0] + ', 0)')
    .call(y_ax);
  const x_called = svg_obj.append('g')
    .attr('class', 'axis')
    .attr('transform', 'translate(0,' + yscale.range()[0] + ')')
    .call(x_ax);

  if (yscale.domain().length > 2) update_axis(y_called);
  
  svg_obj.append('text')
    .attr('class', 'axis_label')
    .attr('text-anchor', 'middle')
    .attr('x', xmid)
    .attr('y', yscale.range()[0]+40)
    .text(xlabel)
    .attr('fill', base_color);

  svg_obj.append('text')
    .attr('class', 'axis_label')
    .attr('text-anchor', 'middle')
    .attr('x', xscale.range()[0])
    .attr('y', yscale.range()[1]-10)
    .text(ylabel)
    .attr('fill', base_color);

  svg_obj.append('text')
    .attr('class', 'axis_label')
    .attr('text-anchor', 'middle')
    .attr('x', xmid)
    .attr('y', yscale.range()[1]-20)
    .text(title)
    .attr('fill', base_color);
}

// GENOME BROWSER //

function gff_parse(r) {
  return {
    'chromosome': r[0], 
    'type': r[2], 
    'start': parseInt(r[3]), 
    'end': parseInt(r[4]),
    'strand': r[6],
    'phase': r[7],
    'attributes': r[8]
  }
}

// https://gist.github.com/tophtucker/62f93a4658387bb61e4510c37e2e97cf
function measureText(string, fontSize = 10) {
  const widths = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.2796875,0.2765625,0.3546875,0.5546875,0.5546875,0.8890625,0.665625,0.190625,0.3328125,0.3328125,0.3890625,0.5828125,0.2765625,0.3328125,0.2765625,0.3015625,0.5546875,0.5546875,0.5546875,0.5546875,0.5546875,0.5546875,0.5546875,0.5546875,0.5546875,0.5546875,0.2765625,0.2765625,0.584375,0.5828125,0.584375,0.5546875,1.0140625,0.665625,0.665625,0.721875,0.721875,0.665625,0.609375,0.7765625,0.721875,0.2765625,0.5,0.665625,0.5546875,0.8328125,0.721875,0.7765625,0.665625,0.7765625,0.721875,0.665625,0.609375,0.721875,0.665625,0.94375,0.665625,0.665625,0.609375,0.2765625,0.3546875,0.2765625,0.4765625,0.5546875,0.3328125,0.5546875,0.5546875,0.5,0.5546875,0.5546875,0.2765625,0.5546875,0.5546875,0.221875,0.240625,0.5,0.221875,0.8328125,0.5546875,0.5546875,0.5546875,0.5546875,0.3328125,0.5,0.2765625,0.5546875,0.5,0.721875,0.5,0.5,0.5,0.3546875,0.259375,0.353125,0.5890625]
  const avg = 0.5279276315789471
  return string
    .split('')
    .map(c => c.charCodeAt(0) < widths.length ? widths[c.charCodeAt(0)] : avg)
    .reduce((cur, acc) => acc + cur) * fontSize
}

class SimpleGenomeBrowser {

  constructor(gff_file, starting_domain, svg_obj, w, h) {
    const self = this;
    this.svg = svg_obj;
    this.w = w;
    this.h = h;
    this.domain = starting_domain;
    this.data_loaded = false
    this.x_range = [layout.inner_buf, self.w-layout.inner_buf];
    this.midpoint = (this.x_range[1]+this.x_range[0])/2;
    d3.text(gff_file).then(function(tdata) {
      self.data = d3.tsvParseRows(tdata.split('\n').filter((line) => (!line.startsWith('#'))).join('\n'), gff_parse);
      self.data = self.data.filter((d) => d.type=='CDS');
      for (let row of self.data) {
        let ats = {}
        row.attributes.split(';').forEach(function(pair) {
          let keyVal = pair.split('=');
          ats[keyVal[0]] = keyVal[1];
        })
        row['gene'] = ats['gene'];
        row['LocusID'] = ats['locus_tag'] || 'NA';
        row['product'] = ats['product'] || 'NA';
        row['gene_index'] = locusID_list.indexOf(row['LocusID']);
      }
      self.data_loaded = true;
      self.data_map = make_locusID_to_row_map(self.data);

      self.drag_zone = self.svg.append('rect')
        .attr('width', self.w)
        .attr('height', self.h)
        .attr('fill', background_color);

      self.zoom_in_g = self.svg.append('g')
        .attr('class', 'zoom_thing')
        .on('click', function() {
          self.display_region(self.zoom_in());
        })
      
      self.zoom_out_g = self.svg.append('g')
        .attr('class', 'zoom_thing')
        .on('click', function() {
          self.display_region(self.zoom_out());
        })
      
      self.zoom_in_g.append('circle')
        .attr('cx', self.midpoint+17)
        .attr('cy', 7)
        .attr('r', 6)
        .attr('stroke', base_color)

      self.zoom_in_g.append('line')
        .attr('x1', self.midpoint+17-3)
        .attr('x2', self.midpoint+17+3)
        .attr('y1', 7)
        .attr('y2', 7)
        .attr('stroke', base_color)

      self.zoom_in_g.append('line')
        .attr('x1', self.midpoint+17)
        .attr('x2', self.midpoint+17)
        .attr('y1', 7-3)
        .attr('y2', 7+3)
        .attr('stroke', base_color)

      self.zoom_out_g.append('circle')
        .attr('cx', self.midpoint-17)
        .attr('cy', 7)
        .attr('r', 6)
        .attr('stroke', base_color)

      self.zoom_out_g.append('line')
        .attr('x1', self.midpoint-17-3)
        .attr('x2', self.midpoint-17+3)
        .attr('y1', 7)
        .attr('y2', 7)
        .attr('stroke', base_color)

      self.x_scale = d3.scaleLinear().range(self.x_range).domain(self.domain);
      self.x_axis = d3.axisTop(self.x_scale).ticks(6);
      self.x_ax_element = self.svg.append('g')
        .attr('class', 'axis')
        .attr('transform', 'translate(0, 20)')
        .call(self.x_axis);

      self.gene_g = self.svg.append('g');
 
      self.dragAction = d3.drag()
        .on('start', function(e) {
          self.drag_start_domain = self.x_scale.domain();
          self.drag_start = self.x_scale.invert(e.x);
          self.drag_start_mouse = e.x;
          self.tmp_scale = d3.scaleLinear().range(self.x_range).domain(self.domain);
          self.tmp_axis = d3.axisTop(self.tmp_scale).ticks(6);
          self.x_ax_element.remove();
          self.x_ax_element = self.svg.append('g')
            .attr('class', 'axis')
            .attr('transform', 'translate(0, 33)')
            .call(self.tmp_axis);
        })
        .on('drag', function(e) {
          const x_pos = self.x_scale.invert(e.x);
          const x_change = x_pos-self.drag_start;
          const x_change_mouse = e.x-self.drag_start_mouse;
          self.tmp_scale.domain([self.drag_start_domain[0]-x_change, self.drag_start_domain[1]-x_change])
          self.x_ax_element.call(self.tmp_axis.scale(self.tmp_scale));
          self.gene_g.attr('transform', 'translate('+String(x_change_mouse)+',0)');
        })
        .on('end', function(e) {
          self.domain = self.tmp_scale.domain();
          self.x_scale = self.tmp_scale;
          self.x_axis = self.tmp_axis;
          self.display_region();
        })
      self.drag_zone.call(self.dragAction);
      self.display_region();
    }).catch(function(error) {
      console.log('Error loading gff file', error);
    })
  }

  zoom_in() {
    this.domain_wid = this.domain[1]-this.domain[0];
    return [Math.max(this.domain[0]+this.domain_wid/4, 0), Math.min(this.domain[1]-this.domain_wid/4, genome_len)];
  }

  zoom_out() {
    this.domain_wid = this.domain[1]-this.domain[0];
    return [Math.max(this.domain[0]-this.domain_wid, 0), Math.min(this.domain[1]+this.domain_wid, genome_len)];
  }

  display_region(new_region=null) {
    if (new_region) this.domain = new_region;
    self = this;
    this.gene_g.remove();
    this.gene_g = this.svg.append('g');
    this.x_scale.domain(self.domain);
    this.x_ax_element.attr('transform', 'translate(0, 33)').call(self.x_axis);
    this.gene_g.selectAll('.sbg_gene')
      .data(self.data.filter((d) => (d.start < self.domain[1]) && (d.end > self.domain[0])))
      .enter()
      .append('g')
        .attr('class', 'sgb_gene')
        .attr('opacity', 0.8)
        .on('mouseover', function(e, d) {
          d3.select(this).attr('opacity', 1);
          d3.select(this).selectAll('polygon').attr('stroke', hover_color).attr('stroke-width', 5);
          const gene_index = locusID_list.indexOf(d.LocusID);
          hover_gene(gene_index);
        })
        .on('mouseout', function (e, d) {
          d3.select(this).attr('opacity', 0.8);
          if (d.LocusID && (d.LocusID == locusID_list[row_clicked])) {
            d3.select(this).selectAll('polygon').attr('stroke', click_color).attr('stroke-width', 5);
          } else {
            d3.select(this).selectAll('polygon').attr('stroke', base_color).attr('stroke-width', 1);
          }
          hover_off();
        })
        .on('click', function (e, d) {
          d3.selectAll('polygon').attr('stroke', base_color).attr('stroke-width', 1);
          d3.select(this).selectAll('polygon').attr('stroke', click_color).attr('stroke-width', 5);
          const gene_index = locusID_list.indexOf(d.LocusID);
          click_gene(gene_index, false);
        })
        .html((d) => self.make_gene_display(d, self))

  }


  make_gene_display(d, self) {
    const clicked = (d.LocusID && (d.LocusID == locusID_list[row_clicked]))
    const left = self.x_scale(d.start);
    const right = self.x_scale(d.end);
    const region_bp = self.x_scale.domain()[1]-self.x_scale.domain()[0];
    const width = right-left;
    const height = Math.max(Math.min(30, 1000000/region_bp), 20);
    const halfHeight = height / 2;
    const chevron_size = (width < 10) ? 0 : Math.min(width/4, 20);
    const top = 41
    const points = d.strand === '-' ? `${left},${top+halfHeight} ${left+chevron_size},${top+height} ${left+width},${top+height} ${left+width},${top} ${left+chevron_size},${top}` : `${right},${top+halfHeight} ${right-chevron_size},${top+height} ${right-width},${top+height} ${right-width},${top} ${right-chevron_size},${top}`;

    const fontsizes = [10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30]
    const textBuf = 2.5
    let label = d.gene || d.LocusID
    let fontsize = fontsizes[0]
    let labelsize = measureText(label, fontsize)
    let labelVisible = (labelsize+2*textBuf+chevron_size < right-left)
    if (labelVisible) {
      for (let f of fontsizes) {
        labelsize = measureText(label, f)
        if (labelsize+2*textBuf+chevron_size < right-left) {
          fontsize = f
        } else {
          break
        }
      }
    }
    const x_pos = d.strand === '+' ? left+textBuf : left+textBuf+chevron_size;
    const y_pos = top+height-textBuf-2;
    const stroke = clicked ? click_color : base_color;
    const strokeWid = clicked ? 5 : 1;
    const chev = `<polygon points="${points}" stroke=${stroke} fill="#333" stroke-width=${strokeWid} />`
    const label_use = labelVisible ? `<text x=${x_pos} y=${y_pos} fill="#FFF">${label}</text>` : '';
    // getting tnseq data
    if (d.gene_index>-1) {
      const tnseq_row = corr_data[d.gene_index];
      const anc_color = tnseq_row['Anc_color']; 
      const focal_color = tnseq_row[focal_pop+'_color'];
      const anc_block = `<rect x=${left} width=${right-left} y=${top+height+10}  height=10 fill="${anc_color}" />`;
      const focal_block = `<rect x=${left} width=${right-left} y=${top+height+25}  height=10 fill="${focal_color}" />`;
      return chev+label_use+anc_block+focal_block;
    } else {
      return chev+label_use;
    }
  }
}

// ALLELE DATA STUFF //

const tf = function(d) {
  //tooltip formatter for allele data
  const g = d.Gene || '';
  return '<h2>'+g+'</h2><p>'+String(d.Location)+' '+String(d.Allele)+'</p><p>'+String(d.Type)+'</p>';
}

// DRAWING
//https://stackoverflow.com/questions/7812514/drawing-a-dot-on-html5-canvas//
function drawPixel(canvasData, canvasWidth, canvasHeight, x, y, r, g, b, a) {
  if ((x < canvasWidth) && (x >= 0) && (y < canvasHeight) && (y >= 0)) {
    let index = (x + y * canvasWidth) * 4;

    canvasData.data[index + 0] += r;
    canvasData.data[index + 1] += g;
    canvasData.data[index + 2] += b;
    canvasData.data[index + 3] += a;
  }
}

function clearImageData(canvasData) {
  for (let i=0; i<canvasData.data.length; i++) {
    canvasData.data[i]=255;
  }
}


class WOW_line_plot {
  /* Plots lines */
  constructor(dimensions, data, parent, xscale, yscale, xlabel, ylabel, 
              pixel_hover_range, tooltip_dat, tooltip_formatter, hover_on) {
    const self = this;
    this.data = data;
    this.xscale = xscale;
    this.yscale = yscale;
    this.tooltip_dat = tooltip_dat;
    this.tooltip_formatter = tooltip_formatter;
    this.hover_on = hover_on;
    this.linefunc = d3.line().x(function(d) { return xscale(d.x); }).y(function(d) { return self.yscale(d.y); })
    this.pixel_hover_range = pixel_hover_range;
    [this.left, this.top, this.w, this.h] = dimensions;
    this.drawing_w = Math.abs(this.xscale.range()[1]-this.xscale.range()[0]);
    this.drawing_h = Math.abs(this.yscale.range()[1]-this.yscale.range()[0]);
    this.lines = this.data.map((d) => this.make_line(d));

    this.canvas = parent.append('canvas')
      .attr('class', 'WOW_canvas')
      .style('left', this.left)
      .style('top', this.top)      
      .attr('width', this.w)
      .attr('height', this.h);

    this.svg = parent.append('svg')
      .attr('class', 'WOW_svg')
      .style('left', this.left)
      .style('top', this.top)
      .attr('width', this.w)
      .attr('height', this.h);

    if (this.hover_on) {
      standard_axes(this.svg, this.xscale, this.yscale, xlabel, ylabel);
      this.hoverMap = [];
      this.hoverMapDistances = [];
      for (let i=0; i<this.w; i++) {
        this.hoverMap.push([]);
        this.hoverMapDistances.push([]);
        for (let j=0; j<this.h; j++) {
          this.hoverMap[i].push(null); 
          this.hoverMapDistances[i].push(pixel_hover_range+1);
        }
      }

      this.svg
        .on('mousemove', function(e, d) {
          let [mx, my] = [e.offsetX, e.offsetY];
          //console.log(mx, my);
          if ((mx < self.w-1) && (my < self.h-1)) {
            let hover_el = self.hoverMap[mx+1][my+1];
            // first just plot the line quick
            self.hoverline.attr('d', self.lines[hover_el] || '');
            if (tooltip_dat && hover_el) {
              // then send round trip to corr_mouseover to plot all lines of the same gene
              const gene_index = locusID_list.indexOf(tooltip_dat[hover_el].LocusID);
              hover_gene(gene_index);
              show_tooltip(e.pageX, e.pageY, tooltip_formatter(tooltip_dat[hover_el]));
            } else {
              hide_tooltip();
            }
          }
        })
        .on('mouseleave', function(e, d) {
          hide_tooltip();
          self.hoverline.attr('d', '');
          hover_off();
        })
        .on('click', function(e, d) {
          let [mx, my] = [e.offsetX, e.offsetY];
          if ((mx < self.w-1) && (my < self.h-1)) {
            let hover_el = self.hoverMap[mx+1][my+1];
            // first just plot the line quick
            self.clickline.append('path').attr('d', self.lines[hover_el] || '');
            if (tooltip_dat && hover_el) {
              // then send round trip to corr_mouseover to plot all lines of the same gene
              const gene_index = locusID_list.indexOf(tooltip_dat[hover_el].LocusID);
              click_gene(gene_index);
            }
          }
        });
    }

    this.ctx = this.canvas.node().getContext('2d');
    this.canvasData = this.ctx.getImageData(0, 0, this.w, this.h);
    
    for (let i=0; i<this.data.length; i++) {
      this.draw_one_line(i, this.data[i]);
    }

    this.ctx.putImageData(this.canvasData, 0, 0);

    this.hoverline = this.svg.append('path')
      .attr('class', 'WOW_hover_path');
      
    this.clickline = this.svg.append('g')
      .attr('class', 'WOW_click_path');
  }

  update_hoverMap() {
    for (let i=0; i<this.w; i++) {
      this.hoverMap.push([]);
      this.hoverMapDistances.push([]);
      this.drawnMap.push([]);
      for (let j=0; j<this.h; j++) {
        this.hoverMap[i].push(null); 
        this.hoverMapDistances[i].push(pixel_hover_range);
        this.drawnMap[i].push([]);
      }
    }  
  }

  make_line(datum) {
    let path = 'M'+String(this.xscale(datum[0][0]))+','+String(this.yscale(datum[1][0]));
    for (let i=1; i<datum[0].length; i++) {
      path += 'L'+String(this.xscale(datum[0][i]))+','+String(this.yscale(datum[1][i]));
    }
    return path;
  }

  hover_check_x(row_index, xi, pixel_x, pixel_y) {
    if ((pixel_x+xi < this.w) && (pixel_x+xi >= 0) && (pixel_y < this.h) && (pixel_y >= 0)) {
      if (Math.abs(xi)<this.hoverMapDistances[pixel_x+xi][pixel_y]) {
        this.hoverMapDistances[pixel_x+xi][pixel_y] = Math.abs(xi);
        this.hoverMap[pixel_x+xi][pixel_y] = row_index;
      }
    }
  }

  hover_check_y(row_index, yi, pixel_x, pixel_y) {
    if ((pixel_x < this.w) && (pixel_x >= 0) && (pixel_y+yi < this.h) && (pixel_y+yi >= 0)) {
      if (Math.abs(yi)<this.hoverMapDistances[pixel_x][pixel_y+yi]) {
        this.hoverMapDistances[pixel_x][pixel_y+yi] = Math.abs(yi);
        this.hoverMap[pixel_x][pixel_y+yi] = row_index;
      }
    }
  }

  round_pixel(x, y) {
    return [Math.round(this.xscale(x)), Math.round(this.yscale(y))];
  }

  draw_segment(row_index, x1raw, y1raw, x2raw, y2raw, r=255, g=255, b=255, a=150, no_aliasing=false) {
    // https://en.wikipedia.org/wiki/Digital_differential_analyzer_(graphics_algorithm)
    // could do https://en.wikipedia.org/wiki/Xiaolin_Wu%27s_line_algorithm
    // if aliasing is a problem
    // only works with x2>x1
    if (x2raw < x1raw) {
      let tmp = [x1raw, y1raw];
      [x1raw, y1raw] = [x2raw, y2raw];
      [x2raw, y2raw] = tmp;
    }
    const [x1, y1] = this.round_pixel(x1raw, y1raw);
    const [x2, y2] = this.round_pixel(x2raw, y2raw);
    const xd = x2-x1;
    const yd = y2-y1;
    const m = yd/xd; //slope
    let step;
    const is_steep = Math.abs(m) > 1;
    if (is_steep) { // primarily moving in y
      step = Math.abs(yd);
    } else {
      step = Math.abs(xd);
    }
    //console.log(x1, y1, x2, y2, step, xd, yd);
    const x_step = xd/step;
    const y_step = yd/step;
    let x = x1;
    let y = y1;
    let pixel_x, pixel_y, pixel_x_low, pixel_x_high, pixel_y_low, pixel_y_high;
    if (is_steep) {
      for (let i=0; i<=step; i++) {
        pixel_x = Math.round(x);
        pixel_y = y;
        if (no_aliasing) {
          drawPixel(this.canvasData, this.w, this.h, pixel_x, pixel_y, r, g, b, a);
        } else {
          // rough anti-aliasing, but weird with overlapping lines
          pixel_x_low = Math.floor(x);
          pixel_x_high = Math.ceil(x);
          drawPixel(this.canvasData, this.w, this.h, pixel_x_low, pixel_y, r, g, b, a*(1-(x-pixel_x_low)));
          drawPixel(this.canvasData, this.w, this.h, pixel_x_high, pixel_y, r, g, b, a*(1-(pixel_x_high-x)));
        }
        if (this.hover_on) {
        // Allow hovering up to pixel_hover_range pixels away from the primary axis (rough approach)
          for (let xi=-1*this.pixel_hover_range; xi<=this.pixel_hover_range; xi++) {
            this.hover_check_x(row_index, xi, pixel_x, pixel_y)
          }
        }
        x += x_step;
        y += y_step;
      }
    } else {
      for (let i=0; i<=step; i++) {
        pixel_x = x;
        pixel_y = Math.round(y);
        if (no_aliasing) {
          drawPixel(this.canvasData, this.w, this.h, pixel_x, pixel_y, r, g, b, a);
        } else {
          // rough anti-aliasing, but weird with overlapping lines
          pixel_y_low = Math.floor(y);
          pixel_y_high = Math.ceil(y);
          drawPixel(this.canvasData, this.w, this.h, pixel_x, pixel_y_low, r, g, b, a*(1-(y-pixel_y_low)));
          drawPixel(this.canvasData, this.w, this.h, pixel_x, pixel_y_high, r, g, b, a*(1-(pixel_y_high-y)));
        }
        if (this.hover_on) {
          // Allow hovering up to pixel_hover_range pixels away from the primary axis (rough approach)
          for (let yi=-1*this.pixel_hover_range; yi<=this.pixel_hover_range; yi++) {
            this.hover_check_y(row_index, yi, pixel_x, pixel_y);
          }
        }
        x += x_step;
        y += y_step;
      }
    }
  }

  draw_one_line(row_index, line) {
    for (let i=0; i<line[0].length-1; i++) {
      this.draw_segment(row_index, line[0][i], line[1][i], line[0][i+1], line[1][i+1]);
    }
  }
}

function one_allele_graph(population, left, top, wid, high, focal) {
  const short_pop = all_pops[formatted_pops.indexOf(population)];
  const allele_div = d3.select('#data_div').append('div')
    .attr('class', 'graph_holder')
    .style('left', left)
    .style('top', top)
    .style('width', wid)
    .style('height', high)

  if (focal) allele_div.attr('id', 'focal_graph');
  let wowplot;
  if (focal) {
    wowplot = new WOW_line_plot([0, 0, wid, high], freq_data[short_pop], allele_div, gens_x, freq_y, 'Generations', 'Allele Freq.', pixel_hover_range=5, tooltip_dat=allele_data[short_pop], tooltip_formatter=tf, hover_on=true);
  } else {
    const overlay = allele_div.append('div')
      .data([population])
      .attr('class', 'graph_overlay')
      .attr('left', 0)
      .attr('top', 0)
      .style('width', wid)
      .style('height', high)
      .on('click', function() { new_focal(population); })
      .on('mouseover', function(e, d) {
        overlay_tooltip(e.pageX, e.pageY, formatted_pops.indexOf(d));
        d3.selectAll('.s_traj').classed('s_traj_hovered', (td) => td.Population==d)
      })
      .on('mouseout', function() {
        hide_tooltip();
        d3.selectAll('.s_traj').classed('s_traj_hovered', false)
      });
    pop_overlays.push(overlay.append('p').html(population));
    const tmp_freq_y = d3.scaleLinear().range([layout.small_h-layout.small_buf, layout.small_buf]).domain([0, 1]);
    wowplot = new WOW_line_plot([0, 0, wid, high], freq_data[short_pop], allele_div, right_halfgens_x, tmp_freq_y, 'Generations', 'Allele Freq.', pixel_hover_range=1, tooltip_dat=allele_data[short_pop], tooltip_formatter=tf, hover_on=false);
  }
  
  wowplot.gene_index_to_muts = allele_data_map[short_pop];
  allele_graphs.push(wowplot);
}

function make_allele_graphs() {
  one_allele_graph(focal_pop, layout.focal_allele_left, layout.focal_allele_top, 
                   layout.focal_allele_width, layout.focal_allele_height, true);
  for (let i=1; i<all_pops.length; i++) {
    one_allele_graph(formatted_pops[i], 
                     layout.small_allele_left+Math.floor((i-1)/6)*layout.small_allele_width/2, 
                     layout.small_allele_top+((i-1)%6)*layout.small_h, 
                     layout.small_allele_width/2, 
                     layout.small_h, false);
  }
  d3.select('#data_div').style('display', 'block');
  d3.select('#loading_text').remove();
}

function new_focal(pop) {
  focal_pop = pop;
  d3.select('#focal_graph').remove();
  one_allele_graph(pop, layout.focal_allele_left, layout.focal_allele_top, 
                   layout.focal_allele_width, layout.focal_allele_height, true);
  d3.select('#corr_graphs_div').remove();
  make_corr_graphs();
  d3.selectAll('.s_traj')
    .attr('stroke', (d) => d.Population == focal_pop ? click_color : base_color)
    .attr('stroke-width', (d) => d.Population == focal_pop ? 3 : 1)
    .filter((d) => d.Population == focal_pop).raise();
}

// FITNESS GRAPH STUFF //

function make_fitness_graph(fit_data) {

  //reformatting data
  const gens = fit_data.map((o) => parseInt(o.Generation))
  let data_reformat = []
  for (let p of Object.keys(fit_data[0])) {
    if (p != 'Generation') data_reformat.push({'Population': format_name(p), 'Fitness': fit_data.map((o) => parseFloat(o[p]))})
  }

  const fit_div = d3.select('#data_div').append('div')
    .attr('class', 'graph_holder')
    .style('left', layout.fit_div_left)
    .style('top', layout.fit_div_top)
    .style('width', layout.fit_div_width)
    .style('height', layout.fit_div_height);

  
  const fitness_svg_obj = fit_div.append('svg')
    .attr('class', 'WOW_svg')
    .style('left', 0)
    .style('top', 0)
    .style('width', layout.fit_div_width)
    .attr('height', layout.fit_div_height);
    
  standard_axes(fitness_svg_obj, right_gens_x, fit_y, 'Generations', 'Fitness');

  fitness_svg_obj.selectAll('.s_traj')
    .data(data_reformat)
    .enter()
    .append('path')
      .attr('class', 's_traj')
      .attr('fill', 'none')
      .attr('stroke', (d) => d.Population == focal_pop ? click_color : base_color)
      .attr('stroke-width', (d) => d.Population == focal_pop ? 3 : 1)
      .attr('d', function(d) {
        var traj_d = [];
        let i = 0;
        for (let gen of gens) {
          traj_d.push({'x': gen, 'y': d['Fitness'][i]})
          i += 1
        }
        //console.log(traj_d);
        return s_traj(traj_d); 
      })
      .on('mouseover', function(e, d) {
        show_tooltip(e.pageX, e.pageY, '<h2>'+d.Population+' </h2>');
        d3.selectAll('.s_traj').filter((d) => d.Population == focal_pop).raise();
        d3.select(this).raise()
        d3.selectAll('.graph_overlay').classed('go_hovered', (td) => td==d.Population);
      })
      .on('mouseout', function(e, d) {
        hide_tooltip();
        d3.selectAll('.graph_overlay').classed('go_hovered', false);
        d3.selectAll('.s_traj').filter((d) => d.Population == focal_pop).raise();
      })
      .on('click', function(e, d) { 
        new_focal(d.Population);
      })
      .filter((d) => d.Population == focal_pop).raise();
}

// TNSEQ DATA GRAPHS //

function overlay_tooltip(x, y, i) {
  if (row_clicked) {
    let p = '';
    for (let mutation_row of allele_graphs[i].gene_index_to_muts[row_clicked]) {
      p += allele_graphs[i].tooltip_formatter(allele_graphs[i].tooltip_dat[mutation_row]);
    }
    if (p!='') show_tooltip(x, y, p);
  }
}

// interaction with allele graphs
function set_allele_lines(which_line, row_index, clear=false) {
  function one_allele_lines(a, which_line, row_index, clear) {
    if (which_line == 'clickline') a[which_line].selectAll('path').remove();
    a.hoverline.attr('d', '');
    if ((~clear) && (row_index) && (row_index != -1)) {
      a[which_line].selectAll('path')
      let path = ''
      for (let mutation_row of a.gene_index_to_muts[row_index]) {
        if (which_line == 'clickline') {
          a[which_line].append('path')
            .attr('d', a.lines[mutation_row])
            .on('mousemove', function(e) {
                e.stopPropagation();
                a.hoverline.attr('d', '');
                show_tooltip(e.pageX, e.pageY, a.tooltip_formatter(a.tooltip_dat[mutation_row]));
              })
            .on('mouseleave', hide_tooltip);
        } else {
          path += a.lines[mutation_row];
        }
      }
      if (which_line == 'hoverline') a[which_line].attr('d', path);
    }
  }
  for (let a of allele_graphs) {
    one_allele_lines(a, which_line, row_index, clear);
  }
  one_allele_lines(spagplot, which_line, row_index, clear);
}

function update_axis(ax) {
  ax.selectAll('text').attr('fill', (d) => d==focal_pop ? hover_color : base_color)
}

function one_corr_graph(svg_obj, xcol, ycol, xscale, yscale) {
  svg_obj.selectAll('.corr_group').append('circle')
    .attr('class', 'corr_point')
    .attr('r', 2)
    .attr('fill', base_color)
    .attr('cx', (d) => xscale(d[xcol]))
    .attr('cy', (d) => yscale(d[ycol]));
  
}

function make_tn_line(datum) {
  if (!datum) return '';
  let p = '';
  const reps = [['Anc', 'rep1'], ['Anc', 'rep2'], [focal_pop, 'rep1'], [focal_pop, 'rep2']];
  for (let i=0; i<4; i++) {
    let r = reps[i]
    p += 'M'+String(tn_x_scales[i](0))+','+String(log_freq_y(datum[r[0]+'_T0_'+r[1]]));
    for (let t=1; t<5; t++) {
      p += 'L'+String(tn_x_scales[i](t))+','+String(log_freq_y(datum[r[0]+'_T'+String(t)+'_'+r[1]]));
    }
  }
  return p;
}

function pop_points_on(hover_row_index=null, turn_on_hover_points=null) {
  // hovering
  d3.selectAll('.overlay_point_hover').attr('opacity', 0);
  if (hover_row_index && (hover_row_index != -1)) {
    d3.selectAll('.pop_point_hover_rep1')
      .attr('opacity', 1)
      .attr('cx', (p) => corr_x3(corr_data[hover_row_index][p+'_rep1']));
    d3.selectAll('.pop_point_hover_rep2')
      .attr('opacity', 1)
      .attr('cx', (p) => corr_x3(corr_data[hover_row_index][p+'_rep2']));
    if (turn_on_hover_points) {
      d3.select('#corr_point_hover_anc')
        .attr('opacity', 1)
        .attr('cx', corr_x1(corr_data[hover_row_index]['Anc_rep1']))
        .attr('cy', corr_y(corr_data[hover_row_index]['Anc_rep2']))
        .attr('locusID', locusID_list[hover_row_index]);

      d3.select('#corr_point_hover_focal')
        .attr('opacity', 1)
        .attr('cx', corr_x2(corr_data[hover_row_index][focal_pop+'_rep1']))
        .attr('cy', corr_y(corr_data[hover_row_index][focal_pop+'_rep2']))
        .attr('locusID', locusID_list[hover_row_index]);

      d3.select('#corr_point_hover_comp')
        .attr('opacity', 1)
        .attr('cx', corr_x3(corr_data[hover_row_index]['Anc_ave']))
        .attr('cy', corr_y(corr_data[hover_row_index][focal_pop+'_ave']))
        .attr('locusID', locusID_list[hover_row_index]);
    }

  
    d3.select('#tn_line_hover_1')
      .attr('opacity', 1)
      .attr('d', make_tn_line(count_data[count_data_map[hover_row_index][0]]));
    d3.select('#tn_line_hover_2')
      .attr('opacity', 1)
      .attr('d', make_tn_line(count_data[count_data_map[hover_row_index][1]]));

  } else if (row_clicked) {
    d3.selectAll('.pop_point_click_rep1')
      .attr('opacity', 1)
      .attr('cx', (p) => corr_x3(corr_data[row_clicked][p+'_rep1']));

    d3.selectAll('.pop_point_click_rep2')
      .attr('opacity', 1)
      .attr('cx', (p) => corr_x3(corr_data[row_clicked][p+'_rep2']));

    d3.select('#corr_point_click_anc')
      .attr('opacity', 1)
      .attr('cx', corr_x1(corr_data[row_clicked]['Anc_rep1']))
      .attr('cy', corr_y(corr_data[row_clicked]['Anc_rep2']));

    d3.select('#corr_point_click_focal')
      .attr('opacity', 1)
      .attr('cx', corr_x2(corr_data[row_clicked][focal_pop+'_rep1']))
      .attr('cy', corr_y(corr_data[row_clicked][focal_pop+'_rep2']));
    
    d3.select('#corr_point_click_comp')
      .attr('opacity', 1)
      .attr('cx', corr_x3(corr_data[row_clicked]['Anc_ave']))
      .attr('cy', corr_y(corr_data[row_clicked][focal_pop+'_ave']));

    d3.select('#tn_line_click_1')
      .attr('opacity', 1)
      .attr('d', make_tn_line(count_data[count_data_map[row_clicked][0]]));
    d3.select('#tn_line_click_2')
      .attr('opacity', 1)
      .attr('d', make_tn_line(count_data[count_data_map[row_clicked][1]]));
  }
}

function hover_points_off() {
  d3.selectAll('.overlay_point_hover').attr('opacity', 0);
  d3.selectAll('.tn_line_hovered').attr('opacity', 0);
}

function click_points_off() {
  d3.selectAll('.overlay_point_click').attr('opacity', 0);
  d3.selectAll('.tn_line_clicked').attr('opacity', 0);
}

function hover_gene(hovered_ind, turn_on_hover_points=true) {
  pop_points_on(hovered_ind, turn_on_hover_points);
  set_allele_lines('hoverline', hovered_ind);
}

function hover_off() {
  hover_points_off();
  set_allele_lines('hoverline', null, true);
}

function click_gene(clicked_ind, update_sgb_display=true) {
  if ((clicked_ind != row_clicked) && (clicked_ind != -1)) { // if it's a new click, style it
    gene_descrip_text.html(gene_list[clicked_ind]+', <i>'+product_list[clicked_ind]+'</i>');
    row_clicked = clicked_ind;
    pop_points_on();
    set_allele_lines('clickline', clicked_ind);
    if (update_sgb_display) {
      const gff_info = sgb.data_map[clicked_ind].map((i) => sgb.data[i]);
      if (gff_info.length == 1) {
        const left = gff_info[0].start;
        const right = gff_info[0].end;
        const wid = right-left;
        sgb.display_region([left-wid*2, right+wid*2]);
      }
    }
  } else {
    gene_descrip_text.html('');
    click_points_off();
    set_allele_lines('clickline', null, true);
    d3.selectAll('polygon').attr('stroke', base_color).attr('stroke-width', 1);
    row_clicked = null;
  }
}

function make_corr_graphs() {
  // adding an index column
  for (let i=0; i< corr_data.length; i++) {
    corr_data[i]['WOW_index'] = i;
  }

  const corr_div = d3.select('#data_div').append('div')
    .attr('class', 'graph_holder')
    .attr('id', 'corr_graphs_div')
    .style('left', layout.corr_div_left)
    .style('top', layout.corr_div_top)
    .style('width', layout.corr_div_width)
    .style('height', layout.corr_div_height);


  const corr_svg_obj = corr_div.append('svg')
    .attr('class', 'WOW_svg')
    .style('left', 0)
    .style('top', 0)
    .attr('width', layout.corr_div_width)
    .attr('height', layout.corr_div_height);

  const corr_svg_overlay = corr_svg_obj.append('g');


  
  corr_svg_overlay.selectAll('.pop_point_hover_rep1')
    .data(formatted_pops)
    .enter()
    .append('circle')
      .attr('class', 'overlay_point_hover pop_point_hover_rep1')
      .attr('r', 4)
      .attr('fill', hover_color)
      .attr('opacity', 0)
      .attr('cy', (p) => pop_y(p)+pop_y_offset);

  corr_svg_overlay.selectAll('.pop_point_hover_rep2')
    .data(formatted_pops)
    .enter()
    .append('circle')
      .attr('class', 'overlay_point_hover pop_point_hover_rep2')
      .attr('r', 4)
      .attr('fill', hover_color)
      .attr('opacity', 0)
      .attr('cy', (p) => pop_y(p)+pop_y_offset);

  corr_svg_overlay.selectAll('.pop_point_click_rep1')
    .data(formatted_pops)
    .enter()
    .append('circle')
      .attr('class', 'overlay_point_click pop_point_click_rep1')
      .attr('r', 4)
      .attr('fill', click_color)
      .attr('opacity', 0)
      .attr('cy', (p) => pop_y(p)+pop_y_offset);

  corr_svg_overlay.selectAll('.pop_point_click_rep2')
    .data(formatted_pops)
    .enter()
    .append('circle')
      .attr('class', 'overlay_point_click pop_point_click_rep2')
      .attr('r', 4)
      .attr('fill', click_color)
      .attr('opacity', 0)
      .attr('cy', (p) => pop_y(p)+pop_y_offset);

  corr_svg_overlay.append('circle')
    .attr('class', 'overlay_point_hover')
    .attr('id', 'corr_point_hover_anc')
    .attr('r', 4)
    .attr('fill', hover_color)
    .attr('opacity', 0);
  corr_svg_overlay.append('circle')
    .attr('class', 'overlay_point_hover')
    .attr('id', 'corr_point_hover_focal')
    .attr('r', 4)
    .attr('fill', hover_color)
    .attr('opacity', 0);
  corr_svg_overlay.append('circle')
    .attr('class', 'overlay_point_hover')
    .attr('id', 'corr_point_hover_comp')
    .attr('r', 4)
    .attr('fill', hover_color)
    .attr('opacity', 0);
  corr_svg_overlay.append('circle')
    .attr('class', 'overlay_point_click')
    .attr('id', 'corr_point_click_anc')
    .attr('r', 4)
    .attr('fill', click_color)
    .attr('opacity', 0)
    .on('click', () => click_gene(row_clicked)); // clears the click
  corr_svg_overlay.append('circle')
    .attr('class', 'overlay_point_click')
    .attr('id', 'corr_point_click_focal')
    .attr('r', 4)
    .attr('fill', click_color)
    .attr('opacity', 0)
    .on('click', () => click_gene(row_clicked)); // clears the click
  corr_svg_overlay.append('circle')
    .attr('class', 'overlay_point_click')
    .attr('id', 'corr_point_click_comp')
    .attr('r', 4)
    .attr('fill', click_color)
    .attr('opacity', 0)
    .on('click', () => click_gene(row_clicked)); // clears the click

  corr_svg_obj.selectAll('.corr_group')
    .data(corr_data)
    .enter()
    .append('g')
      .attr('class', 'corr_group')
      .on('mouseover', function(e, d) {
        d3.select(this).raise().selectAll('.corr_point')
          .attr('fill', hover_color)
          .attr('r', 4);
        hover_gene(locusID_list.indexOf(d.LocusID), false);
        corr_svg_overlay.raise();
        show_tooltip(e.pageX, e.pageY, '<h2>'+d['Gene Name']+' </h2><p>'+d['Protein Product']+'</p>');
      })
      .on('mouseout', function(e, d) {
        d3.select(this).selectAll('.corr_point')
          .attr('fill', base_color)
          .attr('r', 2);
        hover_off();
        hide_tooltip();
      })
      .on('click', function(e, d) {
        click_gene(locusID_list.indexOf(d.LocusID));
      });

  standard_axes(corr_svg_obj, corr_x1, corr_y, 's rep 1', 's rep 2', 'Ancestor TnSeq');
  standard_axes(corr_svg_obj, corr_x2, corr_y, 's rep 1', 's rep 2', focal_pop+' TnSeq');
  standard_axes(corr_svg_obj, corr_x3, corr_y, 's Anc', 's '+focal_pop, '');
    
  one_corr_graph(corr_svg_obj, 'Anc_rep1', 'Anc_rep2', corr_x1, corr_y);
  one_corr_graph(corr_svg_obj, focal_pop+'_rep1', focal_pop+'_rep2', corr_x2, corr_y);
  one_corr_graph(corr_svg_obj, 'Anc_ave', focal_pop+'_ave', corr_x3, corr_y);

  corr_svg_overlay.raise();

  for (let i=0; i<4; i++) {
    tn_x_scales.push(d3.scaleLinear().range([layout.inner_buf+layout.tn_x*i+Math.floor(i/2)*layout.tn_buf*0.5, layout.inner_buf+layout.tn_x*(i+1)-layout.tn_buf+Math.floor(i/2)*layout.tn_buf*0.5]).domain([0,4]))
    standard_axes(corr_svg_obj, tn_x_scales[i], log_freq_y, 'Gens', i==0 ? 'log freq': '', '');
  }

  corr_svg_obj.selectAll('.tn_line')
    .data(count_data.filter((d) => d.Locus_tag=='intergenic'))
    .enter()
    .append('path')
      .attr('class', 'intergenic_line')
      .attr('d', (d) => make_tn_line(d)); 

  corr_svg_obj.append('path')
      .attr('class', 'tn_line_clicked')
      .attr('id', 'tn_line_click_1')
      .attr('opacity', 0);
  corr_svg_obj.append('path')
    .attr('class', 'tn_line_clicked')
    .attr('id', 'tn_line_click_2')
    .attr('opacity', 0);
  corr_svg_obj.append('path')
    .attr('class', 'tn_line_hovered')
    .attr('id', 'tn_line_hover_1')
    .attr('opacity', 0);
  corr_svg_obj.append('path')
    .attr('class', 'tn_line_hovered')
    .attr('id', 'tn_line_hover_2')
    .attr('opacity', 0);

  standard_axes(corr_svg_obj, corr_x3, pop_y, 's', '', '');

  spagplot = new WOW_line_plot([corr_x3.range()[0], pop_y.range()[1], corr_x3.range()[1]-corr_x3.range()[0], pop_y.range()[0]-pop_y.range()[1]], spaghetti_data, corr_div, corr_x3_canvas, pop_y_for_WOW, '', '', pixel_hover_range=1, tooltip_dat=corr_data, tooltip_formatter=((d)=>'<h2>'+d['Gene Name']+' </h2><p>'+d['Protein Product']+'</p>'), hover_on=true);
  spagplot.gene_index_to_muts = spagplot.data.map((row, index) => [index]); // one to one map to LocusID_list by definition
  spagplot.canvas.style('opacity', 0.4);

  //hacky way to make sure clicked gene still shows when we just changed focal population
  if (row_clicked) {
    let hold_click = row_clicked;
    row_clicked = null;
    click_gene(hold_click);
  }
}

// LOADING DATA

function recursive_allele_data_load(pops) {
  if (pops.length > 0) {
    let pop_parse = pops.pop();
    d3.csv('Data_use/LTEE_freqs/'+pop_parse+'_freqs_nononsyn.csv')
      .then(function(raw_data) {
        allele_data[pop_parse] = raw_data;
        let data_reformat = [];
        let gene_index_to_muts = corr_data.map((d) => []);
        let i = 0;
        for (let row of raw_data) {
          let gene_spot = locusID_list.indexOf(row.LocusID);
          if (gene_spot > -1) gene_index_to_muts[gene_spot].push(i);
          data_reformat.push([row['times'].split(';').map((t) => parseInt(t)), row['freqs'].split(';').map((f) => parseFloat(f))]);
          i += 1;
        }
        allele_data_map[pop_parse] = gene_index_to_muts;
        freq_data[pop_parse] = data_reformat;
        recursive_allele_data_load(pops);
      })
      .catch(function(error) {
        console.log(error);
      });  
  } else {
    make_allele_graphs()
  }
}

function load_allele_data() {
  console.log('loading allele data');
  const short_pop = all_pops[formatted_pops.indexOf(focal_pop)];
  let pop_list_copy = all_pops.slice(1, all_pops.length); // excluding Anc
  recursive_allele_data_load(pop_list_copy);
}

function load_fitness_data() {
  d3.csv('Data_use/LTEE_fitness.csv')
    .then(function(raw_data) {
      make_fitness_graph(raw_data);
    })
    .catch(function(error) {
          console.log(error);
    });
}

function make_genome_browser() {
  
  const sgb_div = d3.select('#data_div').append('div')
    .attr('class', 'graph_holder')
    .attr('id', 'sgb_div')
    .style('left', layout.sgb_left)
    .style('top', layout.sgb_top)
    .style('width', layout.sgb_width)
    .style('height', layout.sgb_height);


  const sgb_svg_obj = sgb_div.append('svg')
    .attr('class', 'WOW_svg')
    .style('left', 0)
    .style('top', 0)
    .attr('width', layout.sgb_width)
    .attr('height', layout.sgb_height);

  sgb = new SimpleGenomeBrowser('Data_use/R606.gff3', [4610000, 4620000], sgb_svg_obj, layout.main_col_w+layout.right_col_w, layout.sgb_h);
}

function make_locusID_to_row_map(d) {
  let mapper = corr_data.map((d) => []);
  let i = 0
  for (let row of d) {
    let ind = locusID_list.indexOf(row.LocusID);
    if (ind > -1) {
      mapper[ind].push(i);
    }
    i += 1;
  }
  return mapper;
}

function setup_tooltip() {
  tooltip = d3.select('body').append('div')
    .attr('class', 'WOW_tooltip')
    .html('<h2>yeah</h2><p>uhhuh</p>');
}

function show_tooltip(x, y, text) {
  tooltip.style('left', x-158).style('display', 'block').html(text);
  tooltip.style('top', y-8-tooltip.node().offsetHeight);
}

function hide_tooltip() {
  tooltip.style('display', 'none');
}

// SEARCH BAR //

function make_search_bar() {
  search_div = d3.select('#main_div').append('div')
    .attr('id', 'search_div');
  gene_search = search_div.append('input')
    .attr('id', "gene_searchbar")
    .attr('type', 'search')
    .attr('spellcheck', 'false')
    .attr('autocomplete', 'off');

  const gene_descrip = search_div.append('div')
    .attr('id', "gene_descrip");

  gene_descrip_text = gene_descrip.append('div')
    .attr('id', "gene_descrip_text")
    .html('');

  //type="search" dir="ltr" spellcheck=false autocorrect="off" autocomplete="off" autocapitalize="off"
  const autoCompleteJS = new autoComplete({
      placeHolder: "Search for gene names...",
      selector: '#gene_searchbar',
      data: {src: gene_list},
      events: {
          input: {
              selection: (event) => {
                  const selection = event.detail.selection.value;
                  autoCompleteJS.input.value = '';
                  click_gene(gene_list.indexOf(selection));
                  document.getElementById('gene_searchbar').blur(); //removes focus so the cursor leaves
              }
          }
      }
  });

}

function load_data() {
  console.log('loading tnfit data');
  d3.csv('Data_use/Limdi_fitness_data.csv')
    .then(function(raw_data) {
      corr_data = raw_data; // using a global dataframe here
      for (let row of corr_data) {
        if (row['LocusID']=='') {
          row['LocusID'] = row['Gene Name'];
        }
        for (let p of formatted_pops) {
          //row[p+'_ave'] = (parseFloat(row[p+'_rep1'])+parseFloat(row[p+'_rep2']))/2
          row[p+'_color'] = color_scale(row[p+'_ave']);
        }
        spaghetti_data.push([formatted_pops.map((p) => row[p+'_ave']), pop_y_array]);
      }
      locusID_list = corr_data.map((d) => d['LocusID']);
      gene_list = corr_data.map((d) => d['Gene Name']);
      product_list = corr_data.map((d) => d['Protein Product']);
      // load count data
      d3.csv('Data_use/Limdi_count_data_pooled.csv')
      .then(function(raw_data) {
        count_data = raw_data;
        count_data_map = make_locusID_to_row_map(count_data);
        make_corr_graphs();
        load_allele_data();
        load_fitness_data();
        make_genome_browser();
        setup_tooltip();
        make_search_bar();
      }).catch(function(error) {
        console.log(error);
      });
    })
    .catch(function(error) {
      console.log(error);
    });
}