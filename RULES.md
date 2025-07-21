# Development Rules

## Code Style
- Keep code as simple as possible without compromising functionality
- No comments in code - code should be self-explanatory
- Prefer clear variable names and function names over comments
- Remove all debugging console.log statements in production code
- Keep functions small and focused

## UI Design
- "Stupid simple" UI with no extraneous lines or outlining
- Icon-first mobile UI - prioritize icons over text on small screens
- Minimal visual clutter
- Clean, modern interface

## Development Workflow
- Always commit and push after acceptance
- Use developer account system instead of special dev modes
- Test functionality thoroughly before committing
- Keep changes focused and atomic

## Molecular Visualization
- Use ONLY sphere representation for molecules with van der Waals radii at 0.8 scale
- NO ball-and-stick models
- Display molecules in columns with descriptive names and individual close buttons

## Error Handling
- Provide clear, actionable error messages
- Graceful degradation when features fail
- User-friendly feedback for all interactions

## Performance
- Optimize for mobile devices
- Minimize API calls
- Efficient image processing and analysis 