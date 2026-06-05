import { defineConfig } from 'vitepress'
import { tabsMarkdownPlugin } from 'vitepress-plugin-tabs'
import mathjax3 from "markdown-it-mathjax3";
import footnote from "markdown-it-footnote";
import path from 'path'
import type { ShikiTransformer } from 'shiki'

function getBaseRepository(base: string): string {
  if (!base || base === '/') return '/';
  const parts = base.split('/').filter(Boolean);
  return parts.length > 0 ? `/${parts[0]}/` : '/';
}

const baseTemp = {
  base: 'REPLACE_ME_DOCUMENTER_VITEPRESS',
}

const navTemp = {
  nav: 'REPLACE_ME_DOCUMENTER_VITEPRESS',
}

const sidebarTemp = {
  sidebar: 'REPLACE_ME_DOCUMENTER_VITEPRESS',
}

function collapseAllSidebarGroups(sidebar: any): any {
  const collapseInArray = (arr: any[]) =>
    arr.map((entry) => {
      if (!entry || typeof entry !== "object") return entry
      if (Array.isArray((entry as any).items)) {
        return {
          ...entry,
          collapsible: true,
          collapsed: true,
          items: collapseInArray((entry as any).items),
        }
      }
      return entry
    })
  if (Array.isArray(sidebar)) return collapseInArray(sidebar)
  if (sidebar && typeof sidebar === "object") {
    const out: Record<string, any> = {}
    for (const key of Object.keys(sidebar)) {
      const v = (sidebar as any)[key]
      out[key] = Array.isArray(v) ? collapseInArray(v) : v
    }
    return out
  }
  return sidebar
}

const nav = [
  ...navTemp.nav,
  { component: 'VersionPicker' }
]

// ============================================================================
// PyData a11y-high-contrast syntax highlighting
//
// PyData Sphinx theme uses a11y-high-contrast-light/dark as its default Pygments
// styles. We reproduce those exact token colors via a Shiki transformer that
// remaps github-light/dark output colors.
//
// The key challenge: github-light uses one color (#D73A49) for both keywords
// (using, function, end) and operators (=, +, ==, ||). Pygments distinguishes
// them — keywords are purple (#6730c5), operators are green (#00622f). We split
// them by checking the token text against a known keyword list.
//
// Color mapping (github-light -> a11y-hc-light):
//   #D73A49  keyword/operator  -> #6730c5 (keyword) / #00622f (operator)
//   #005CC5  builtin/number    -> #7f4707 amber
//   #6F42C1  function name     -> #005b82 teal-blue
//   #6A737D  comment           -> #515151 gray
//   #032F62  string            -> #00622f dark green
//   #24292E  default text      -> #080808 near-black
//   #fff     background        -> #fefefe
//
// Color mapping (github-dark -> a11y-hc-dark):
//   #F97583  keyword/operator  -> #dcc6e0 (keyword) / #abe338 (operator)
//   #79B8FF  builtin/number    -> #ffd900 yellow
//   #B392F0  function name     -> #00e0e0 cyan
//   #6A737D  comment           -> #ffd900 yellow
//   #9ECBFF  string            -> #abe338 bright green
//   #E1E4E8  default text      -> #f8f8f2 near-white
//   #24292e  background        -> #2b2b2b
// ============================================================================

const JULIA_KEYWORDS = new Set([
  'using', 'import', 'export', 'module', 'baremodule',
  'function', 'macro', 'return', 'end',
  'if', 'elseif', 'else',
  'for', 'while', 'do', 'in',
  'try', 'catch', 'finally', 'throw',
  'struct', 'mutable', 'abstract', 'type', 'primitive',
  'const', 'local', 'global',
  'let', 'begin', 'quote',
  'where', 'new', 'break', 'continue',
])

const LIGHT_SPAN_MAP: Record<string, string> = {
  '#D73A49': '__SPLIT_LIGHT__', '#d73a49': '__SPLIT_LIGHT__',
  '#005CC5': '#7f4707',  '#005cc5': '#7f4707',
  '#6F42C1': '#005b82',  '#6f42c1': '#005b82',
  '#6A737D': '#515151',  '#6a737d': '#515151',
  '#032F62': '#00622f',  '#032f62': '#00622f',
  '#24292E': '#080808',  '#24292e': '#080808',
}

const DARK_SPAN_MAP: Record<string, string> = {
  '#F97583': '__SPLIT_DARK__', '#f97583': '__SPLIT_DARK__',
  '#79B8FF': '#ffd900',  '#79b8ff': '#ffd900',
  '#B392F0': '#00e0e0',  '#b392f0': '#00e0e0',
  '#6A737D': '#ffd900',  '#6a737d': '#ffd900',
  '#9ECBFF': '#abe338',  '#9ecbff': '#abe338',
  '#E1E4E8': '#f8f8f2',  '#e1e4e8': '#f8f8f2',
}

function applyColorMap(style: string, map: Record<string, string>, text: string): string {
  let s = style
  for (const [from, to] of Object.entries(map)) {
    if (to === '__SPLIT_LIGHT__') {
      s = s.replaceAll(from, JULIA_KEYWORDS.has(text) ? '#6730c5' : '#00622f')
    } else if (to === '__SPLIT_DARK__') {
      s = s.replaceAll(from, JULIA_KEYWORDS.has(text) ? '#dcc6e0' : '#abe338')
    } else {
      s = s.replaceAll(from, to)
    }
  }
  return s
}

const pydataA11yTransformer: ShikiTransformer = {
  name: 'pydata-a11y',
  span(node) {
    const style = node.properties?.style
    if (typeof style !== 'string') return
    const text = (node.children as any[])
      ?.filter(c => c.type === 'text')
      .map(c => c.value)
      .join('')
      .trim() ?? ''
    let s = applyColorMap(style, LIGHT_SPAN_MAP, text)
    s = applyColorMap(s, DARK_SPAN_MAP, text)
    if (s !== style) node.properties.style = s
  },
  pre(node) {
    const style = node.properties?.style
    if (typeof style !== 'string') return
    const isLight = /background-color:#fff(?![0-9a-f])/i.test(style)
                 || /background-color:#ffffff/i.test(style)
    const isDark  = /background-color:#24292e/i.test(style)
    let s = style
    if (isLight) {
      s = s.replace(/#fff(?![0-9a-f])/gi, '#fefefe').replace(/#ffffff/gi, '#fefefe')
      s = applyColorMap(s, LIGHT_SPAN_MAP, '')
    } else if (isDark) {
      s = s.replace(/#24292e/gi, '#2b2b2b')
      s = applyColorMap(s, DARK_SPAN_MAP, '')
    }
    if (s !== style) node.properties.style = s
  },
}

export default defineConfig({
  base: 'REPLACE_ME_DOCUMENTER_VITEPRESS',
  title: 'REPLACE_ME_DOCUMENTER_VITEPRESS',
  description: 'REPLACE_ME_DOCUMENTER_VITEPRESS',
  lastUpdated: true,
  cleanUrls: true,
  outDir: 'REPLACE_ME_DOCUMENTER_VITEPRESS',
  head: [
    ['link', { rel: 'icon', href: 'REPLACE_ME_DOCUMENTER_VITEPRESS_FAVICON' }],
    ['script', {src: `${getBaseRepository(baseTemp.base)}versions.js`}],
    ['script', {src: `${baseTemp.base}siteinfo.js`}]
  ],

  vite: {
    define: {
      __DEPLOY_ABSPATH__: JSON.stringify('REPLACE_ME_DOCUMENTER_VITEPRESS_DEPLOY_ABSPATH'),
    },
    resolve: {
      alias: { '@': path.resolve(__dirname, '../components') }
    },
    optimizeDeps: {
      exclude: [
        '@nolebase/vitepress-plugin-enhanced-readabilities/client',
        'vitepress',
        '@nolebase/ui',
      ],
    },
    ssr: {
      noExternal: [
        '@nolebase/vitepress-plugin-enhanced-readabilities',
        '@nolebase/ui',
      ],
    },
  },

  markdown: {
    math: true,
    config(md) {
      md.use(tabsMarkdownPlugin)
      md.use(mathjax3)
      md.use(footnote)
    },
    theme: { light: 'github-light', dark: 'github-dark' },
    codeTransformers: [pydataA11yTransformer],
  },

  themeConfig: {
    outline: 'deep',
    logo: 'REPLACE_ME_DOCUMENTER_VITEPRESS',
    search: {
      provider: 'local',
      options: { detailedView: true }
    },
    nav,
    sidebar: collapseAllSidebarGroups(sidebarTemp.sidebar as any),
    editLink: 'REPLACE_ME_DOCUMENTER_VITEPRESS',
    socialLinks: [
      { icon: 'github', link: 'https://github.com/sebapersson/SBMLImporter.jl' }
    ],
    footer: {
      message: 'Made with <a href="https://luxdl.github.io/DocumenterVitepress.jl/dev/" target="_blank"><strong>DocumenterVitepress.jl</strong></a><br>',
      copyright: `© Copyright ${new Date().getUTCFullYear()}.`
    }
  }
})
